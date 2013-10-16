function varargout = sw_version()
% SW_VERSION gives the current version of the SpinW code
%

% read file header from sw.m file
%fid = fopen([sw_rootdir 'm_files' filesep '@sw' filesep 'sw.m']);
fid = fopen('sw.m');

% first line
fgets(fid);
% skip comments
fLine = strtrim(fgets(fid));
while numel(fLine)>0 && strcmp(fLine(1),'%')
    fLine = strtrim(fgets(fid));
end
% skip empty lines
while numel(fLine)==0 || ~strcmp(fLine(1),'%')
    fLine = strtrim(fgets(fid));
end
% read version information
verLine = {fLine};
while numel(verLine{end})>0 && strcmp(verLine{end}(1),'%')
    verLine{end+1} = strtrim(fgets(fid)); %#ok<*AGROW>
end
verLine = verLine(1:end-1);
fclose(fid);

% read strings enclosed with $ signs
partStr = {};
for ii = 1:numel(verLine)
    verSel = verLine{ii};
    while sum(verSel=='$') > 1
        [~, verSel] = strtok(verSel,'$'); %#ok<*STTOK>
        [partStr{end+1}, verSel] = strtok(verSel,'$');
    end
    
end

nField = numel(partStr);
fieldName = cell(1,nField);
fieldVal = cell(1,nField);

verStruct = struct;

% extract values and save them into a structure
for ii = 1:nField
    [fieldName{ii}, fieldVal{ii}] = strtok(partStr{ii},':');
    fieldVal{ii} = strtrim(fieldVal{ii}(2:end));
    verStruct.(fieldName{ii}) = fieldVal{ii};
end

if nField == 0
    varargout{1} = [];
end

if nargout < 1
    disp([verStruct.Name verStruct.Version ' (rev ' num2str(verStruct.Revision) ')']);
else
    varargout{1} = verStruct;
end

end