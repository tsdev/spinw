function data = sw_readspec(path)
% data = SW_READSPEC(datapath) reads experimental spin wave dispersion data
% from text file, whose location is defined in path.
%
% Format of the input data file:
% Every line of the data file contains information about an energy scan at
% constant Q. Data consists of floating point numbers separated by space
% (first line can be a header line):
%
%   QH QK QL minE maxE I1 E1 w1 I2 E2 w2 ...
%       where:
% QH    H index of the Q point,
% QK    K index of the Q point,
% QL    L index of the Q point,
% minE  lower boundary of the E scan,
% maxE  upper boundary of the E scan,
% In    intensity of the n-th spin wave mode,
% En    center of the n-th spin wave mode,
% wn    width of the n-th spin wave mode.
% The number of modes in a single line of the data file is unlimited,
% however in every line the number of modes have to be the same. Scans with
% less modes should contain in the end zero intensities.
%
% Before any data line a special line can be inserted that gives the
% measured correlation in square brackets, for axample: [Mxx+Myy]. For the
% formatting of this string, see <a href="matlab:doc sw_parstr">sw_parstr</a>.
% If the measured type of correlation is undefined, unpolarised neutron
% scattering intensity is assumed ([Sperp]). When cross sections measured
% in the Blume-Maleev coordinate system (see <a href="matlab:doc sw_egrid">sw_egrid</a>), the normal to the
% scattering plane has to be also defined. This can be given in a second
% pair of square brackes in the xyz coordinate system, for example: [Myy]
% [1 0 0]. If n is undefined, the default value is [0 0 1].
%
%
% Example input data file (polarised scans in the (0KL) plane):
% QH    QK        QL      ENlim1  ENlim2  I1  EN1       W1    I2  EN2       W2
% [Mxx] [1 0 0]
% 0     1        2.9992   0       15      1    3.7128   1.0   1   8.6778    1.0
% 0     1        2.8993   0       15      1    7.0000   1.0   1   11.1249   1.0
% 0     1        2.7993   0       20      1   13.8576   1.0   0   0.0       0.0
% 0     1        2.6994   0       20      1   17.3861   1.0   0   0.0       0.0
% [Myy] [1 0 0]
% 0     1.0000   2.0000   0       25      1   20.2183   1.0   0   0.0       0.0
% 0     1.1000   2.0000   15      30      1   22.7032   1.0   0   0.0       0.0
% 0     1.2000   2.0000   20      35      1   25.1516   1.0   0   0.0       0.0
%
% See also SW.FITSPEC, SW_PARSTR.
%

if nargin == 0
    help sw_readspec;
    return;
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
    error('sw:sw_readspec:DataFormatError','Wrong input dat format!');
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
                error('sw:sw_readspec:DataFormatError','Wrong input dat format!');
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
        data{polIdx}.I     = dTemp(:,6:3:end)';
        data{polIdx}.E     = dTemp(:,7:3:end)';
        data{polIdx}.w     = dTemp(:,8:3:end)';
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