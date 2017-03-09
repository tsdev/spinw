function dat = sw_readtable(dataSource,delimiter,nHeader)
% reads tabular data
%
% dat = SW_READTABLE(dataSource, {delimiter},{nHeader})
%
% Function reads tabular data that has arbitrary header lines denoted with
% # and the last header line is followed by a column name line. The data
% can be arbitrary combination of strings and numbers. A predefined field
% will be added to the imported data 'MODE', this field contains any
% additional string given between data rows, otherwise empty.
%
% Input:
%
% dataSource    Data source, can be file, url or string.
% delimiter     Delimiter of the data, default is whitespace.
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

dataStr = ndbase.source(dataSource);

% if dataStr == -1
%     error('sw_readtable:FileNotFound',['Data file not found: '...
%         regexprep(dataSource,'\' , '\\\') '!']);
% end

% split string into lines
dataStr = regexp(dataStr, ['(?:' sprintf('\n') ')+'], 'split');
% add '\n' back to the end of lines
dataStr = cellfun(@(C)[C char(10)],dataStr,'UniformOutput',false);

% index into the line number, skip header lines
idxStr = nHeader+1;

% read header
while isempty(dataStr{idxStr}) || dataStr{idxStr}(1) == '#'
    idxStr = idxStr+1;
end

% read column names
cTemp = regexp(strtrim(dataStr{idxStr}), ['(?:', delimiter, ')+'], 'split');

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

% next line
idxStr = idxStr+1;

% load data
% remove header lines
dataStr = dataStr(idxStr:end);
% remove empty lines
dataStr = dataStr(cellfun(@(C)~isempty(C),dataStr));
% find mode string positions
isModeStr = cellfun(@(C)C(1)=='#',dataStr);
modeStr = strtrim([{''} dataStr(isModeStr)]);
% fill out the mode string
modeStrIdx = cumsum(isModeStr)+1;
% remove mode strings from data
dataStr    = dataStr(~isModeStr);
modeStrIdx = modeStrIdx(~isModeStr);
modeStr    = modeStr(modeStrIdx);

% read the first line of data to identify column types string/float
firstLine = regexp(strtrim(dataStr{1}),['(?:', delimiter, ')+'],'split');
datFormat = '';
for ii = 1:nCol
    if isnan(str2double(firstLine{ii}))
        % string type
        datFormat = [datFormat '%s ']; %#ok<AGROW>
    else
        datFormat = [datFormat '%f ']; %#ok<AGROW>
    end
end
datFormat = datFormat(1:(end-1));

% read all data
nDat = numel(dataStr);
if delimiter(1) == ' '
    datTemp = textscan([dataStr{:}],datFormat,nDat);
else
    datTemp = textscan([dataStr{:}],datFormat,nDat,'Delimiter',delimiter);
end

% convert the data into struct format
% create an empty structure
sTemp = [cName;repmat({cell(nDat,1)},1,numel(cName))];
dat = struct(sTemp{:});

% fill in all fields
[dat(:).MODE] = modeStr{:};
for ii = 2:numel(cName)
    colSel = ismember(cSel,cName{ii});
    datMat = [datTemp{colSel}];
    if ~iscell(datMat)
        datMat = mat2cell(datMat,ones(1,nDat),sum(colSel));
    end
    [dat(:).(cName{ii})] = datMat{:};
end

end