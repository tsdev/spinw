function dat = sw_readtable(fName,delimiter,nHeader)
% reads tabular data
%
% dat = SW_READTABLE(fName, {delimiter},{nHeader})
%
% Function reads tabular data that has arbitrary header lines denoted with
% # and the last header line is followed by a column name line. The data
% can be arbitrary combination of strings and numbers. A predefined field
% will be added to the imported data 'MODE', this field contains any
% additional string given between data rows, otherwise empty.
%
% Input:
%
% fName     File name string.
% delimiter Delimiter of the data, default is whitespace.
%
% Output:
%
% dat       Structure with field defined in the data.
%
% Example:
%
% The following file (test.dat) is given as input:
% # TEST DATA
% Q(1) Q(2)        Q(3) ENlim(1) ENlim(2) I(1)  EN(1)  s(1) I(2)   EN(2)   s(2)
% # [Mxx] [1 0 0]
% 0     1        2.9992   0       15      1    3.7128   1.0   1   8.6778    1.0
% 0     1        2.8993   0       15      1    7.0000   1.0   1   11.1249   1.0
% 0     1        2.7993   0       20      1   13.8576   1.0   0   0.0       0.0
% 0     1        2.6994   0       20      1   17.3861   1.0   0   0.0       0.0
% # [Myy] [1 0 0]
% 0     1.0000   2.0000   0       25      1   20.2183   1.0   0   0.0       0.0
% 0     1.1000   2.0000   15      30      1   22.7032   1.0   0   0.0       0.0
% 0     1.2000   2.0000   20      35      1   25.1516   1.0   0   0.0       0.0
%
% The command to import the data
% >> dat = sw_readtable('test.dat');
%
% The dat variable will contain the fields MODE, Q, ENlim, I, EN and s and
% it will have 7 entry. The mode variable will contain '[Mxx] [1 0 0]'
% string for the first 4 entry and '[Myy] [1 0 0]' for the last 3 entry.
% For example the field Q has 3 elements per entry, to extract all Q points
% into a matrix use the command:
% >> Q = reshape([dat(:).Q],3,[])';
%

if nargin == 0
    help sw_readtable
    return
end

if nargin == 1
    delimiter = ' ';
end

if nargin < 3
    nHeader = 0;
end

fid = fopen(fName);

if fid == -1
    error('spinw:sw_readtable:FileNotFound',['Data file not found: '...
        regexprep(fName,'\' , '\\\') '!']);
end

% read header lines given by user
if nHeader > 0
    for ii = 1:nHeader
        fgets(fid);
    end
end

% read header
str = fgets(fid);
while isempty(str) || str(1) == '#'
    str = fgets(fid);
end

% read column names
%cTemp = strsplit(strtrim(str),delimiter);
cTemp = regexp(strtrim(str), ['(?:', delimiter, ')+'], 'split');

cName = {};
mIdx  = cell(1,numel(cTemp));

cSel = cell(1,numel(cTemp));

% read indices from column names
for ii = 1:numel(cTemp)
    sTemp1 = strfind(cTemp{ii},'(');
    sTemp2 = strfind(cTemp{ii},')');
    if ~isempty(sTemp1)
        cIdx = cTemp{ii}((sTemp1+1):(sTemp2-1));
    else
        cIdx = '1';
        sTemp1 = numel(cTemp{ii})+1;
    end
    % remove spaces and find indices
    mIdx{ii} = num2cell(sscanf(cIdx(cIdx~=' '),'%d,'));
    
    hIdx = find(strcmp(cName,cTemp{ii}(1:(sTemp1(end)-1))));
    if isempty(hIdx)
        cName{end+1} = cTemp{ii}(1:(sTemp1(end)-1)); %#ok<AGROW>
        cSel(ii) = cName(end);
    else
        cSel(ii) = cName(hIdx(1));
    end
    
end

cName = [{'MODE'} cName];
nCol  = numel(cSel);

% create an empty structure
sTemp = [cName;repmat({{}},1,numel(cName))];
dat = struct(sTemp{:});

modeStr = ''; 
idx = 1;

str = fgets(fid);

% load data
while str ~= -1
    
    str = strtrim(str);
    
    %lTemp = strsplit(str,delimiter);
    lTemp = regexp(str,['(?:', delimiter, ')+'],'split');
    
    if ~isempty(str) && str(1) == '#'
        modeStr = str;
        str = fgets(fid);
        continue
    end
    
    for ii = 1:nCol
        val = sscanf(lTemp{ii},'%f');
        if isempty(val)
            % cannot index strings
            dat(idx).(cSel{ii}) = lTemp{ii};
        else
            dat(idx).(cSel{ii})(mIdx{ii}{:}) = val;
        end
        
    end
    dat(idx).MODE = modeStr;
    idx = idx + 1;
    str = fgets(fid);
end

fclose(fid);

end