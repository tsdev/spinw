function data = sw_readspec(path)
% read spin wave dispersion data from file
% 
% ### Syntax
% 
% `data = sw_readspec(datapath)`
% 
% ### Description
% 
% `data = sw_readspec(datapath)` reads experimental spin wave dispersion
% data from a text file at the given location. The general format of the
% data file is described in [sw_readtable] with `sw_readspec` requires
% predefined header names and only specific tag strings are allowed. The
% following header line is required:
%  
% ```none
% QH QK QL minE maxE I1 E1 s1 I2 E2 s2 ...
% ```
% where:
%  
% * `QH`        $h$ index of the $Q$ point in rlu,
% * `QK`        $k$ index of the $Q$ point in rlu,
% * `QL`        $l$ index of the $Q$ point in rlu,
% * `minE`      lower boundary of the $E$ scan,
% * `maxE`      upper boundary of the $E$ scan,
% * `In`        intensity of the $n$th spin wave mode,
% * `En`        center of the $n$th spin wave mode, has to be in increasing order,
% * `sn`        standard deviation of the corresponding energy
%  
% The number of modes in a single line of the data file is unlimited,
% however in every line the number of modes have to be the same. Rows with
% less modes should contain zero intensities at the position of the missing
% modes.
%
% {{note `sw_readspec` omits modes that have either zero intensity
% or zero energy.}}
%  
% Before any data line a special tag line can be inserted that gives the
% measured correlation in square brackets, for axample: `'[Mxx+Myy]'`. For
% the formatting of this string, see [sw_parstr]. If the measured type of
% correlation is undefined, unpolarised neutron scattering intensity is
% assumed (equivalent to `'Sperp'`). When cross sections measured in the
% Blume-Maleev coordinate system, see [sw_egrid], the normal to the
% scattering plane has to be also defined. This can be given in a second
% pair of square brackes in the $xyz$ coordinate system, for example: `'[Myy]
% [1 0 0]'`. If $n$ is undefined, the default value is `'[0 0 1]'`.
%  
% ### Examples
%
% Example input data file (polarised scans in the $(0,k,l)$ scattering plane):
% 
% ```none
% QH    QK        QL      ENlim1  ENlim2  I1  EN1       s1    I2  EN2       s2
% [Mxx] [1 0 0]
% 0     1        2.9992   0       15      1    3.7128   1.0   1   8.6778    1.0
% 0     1        2.8993   0       15      1    7.0000   1.0   1   11.1249   1.0
% 0     1        2.7993   0       20      1   13.8576   1.0   0   0.0       0.0
% 0     1        2.6994   0       20      1   17.3861   1.0   0   0.0       0.0
% [Myy] [1 0 0]
% 0     1.0000   2.0000   0       25      1   20.2183   1.0   0   0.0       0.0
% 0     1.1000   2.0000   15      30      1   22.7032   1.0   0   0.0       0.0
% 0     1.2000   2.0000   20      35      1   25.1516   1.0   0   0.0       0.0
% ```
% 
% ### See Also
% 
% [sw_egrid] \| [spinw.fitspec]
%
% *[rlu]: reciprocal lattice unit
%

if nargin == 0
    help sw_readspec
    return
end

fid = fopen(path);

temp  = fgets(fid);
dTemp = sscanf(temp,'%f');
if isempty(dTemp) && ~strcmp(temp(1),'[')
    % The first line was a header read the next one
    temp = fgets(fid);
    dTemp = sscanf(temp,'%f');
end

if ~isempty(dTemp)
    modeStr{1} = 'Sperp';
elseif strcmp(temp(1),'[')
    strPosR = strfind(temp,']')-1;
    modeStr{1} = temp(2:strPosR(1)); %#ok<*AGROW>
else
    error('sw_readspec:DataFormatError','Wrong input dat format!');
end

polIdx = 1;
data   = {};

dTemp = dTemp';
while ~feof(fid)
    temp  = fgets(fid);
    tempF = sscanf(temp,'%f');
    if ~isempty(tempF)
        % Read only lines without text
        dTemp(end+1,:) = tempF';
    end
    if (isempty(tempF) && strcmp(temp(1),'[')) || feof(fid)
        % Store the string that determines which cross section is measured.
        % Remove the trailing stuff after ']' symbol.
        strPosR = strfind(temp,']')-1;
        if numel(strPosR>1)
            strPosL = strfind(temp,'[');
            if numel(strPosL)<2
                error('sw_readspec:DataFormatError','Wrong input data format!');
            end
            nStr = temp(strPosL(2)+1:strPosR(2)-1);
            n    = sscanf(nStr,'%f',[1 3]);
        else
            n    = [0 0 1];
        end
        
        data{polIdx}.n     = n;
        data{polIdx}.Q     = dTemp(:,1:3)';
        data{polIdx}.minE  = dTemp(:,4)';
        data{polIdx}.maxE  = dTemp(:,5)';
        data{polIdx}.E     = dTemp(:,7:3:end)';
        data{polIdx}.I     = dTemp(:,6:3:end)';
        data{polIdx}.sigma = dTemp(:,8:3:end)';
        
        % set intensity zero where energy is zero
        % we won't fit zero energy modes
        data{polIdx}.I(data{polIdx}.E==0) = 0;
        
        % sort intensity, put zero intensities to the end
        [data{polIdx}.I,idx] = sort(data{polIdx}.I,1,'descend');
        data{polIdx}.E       = data{polIdx}.E(sub2ind(size(data{polIdx}.E),idx,repmat(1:size(data{polIdx}.E,2),[2 1])));
        data{polIdx}.sigma   = data{polIdx}.sigma(sub2ind(size(data{polIdx}.E),idx,repmat(1:size(data{polIdx}.E,2),[2 1])));
        
        data{polIdx}.nMode = sum(data{polIdx}.I~=0,1);
        data{polIdx}.corr  = sw_parstr(modeStr{polIdx});
        
        if ~feof(fid)
            modeStr{polIdx+1}  = temp(2:strPosR(1));
        end
        dTemp = [];
        polIdx = polIdx + 1;
        
    end
    
end

fclose(fid);

end