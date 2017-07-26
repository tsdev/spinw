function helpStr = help(name)
% provide help strings in the compiled application
%
% CHANGE THIS TO HELP()
%

if ~isdeployed
    helpStr = '';
    return
end

% get the location of the m-files
[appPath, appName] = sw_apppath;
appStr = [ appName fs 'Source'];

% name without .m extension
name1 = strsplit(name,'.');
% find the right file
[~,fName] = system(['find ' appPath filesep appStr filesep ' -type f | grep "/' name1{1} '\.m"']);

% default return value
helpStr = '';

if ~isempty(fName)
    % extract help text
    helpStr = strtrim(strsplit(fileread(strtrim(fName)),'\n'));
    helpStr(cellfun(@(C)numel(C),helpStr,'uniformoutput',true)==0) = [];
    isComment = cellfun(@(C)C(1)=='%',helpStr);
    lIdx1 = find(isComment,1,'first');
    lIdx2 = find(diff(isComment)==-1,1,'first');
    helpStr = cellfun(@(C)C(2:end),helpStr(lIdx1:lIdx2),'uniformoutput',false);
    helpStr = sprintf('%s\n',helpStr{:});
else
    % try directories
    fName = cell(1,3);
    [~,fName{1}] = system(['find ' appPath filesep appStr filesep ' -type d | grep "/' name1{1} '"']);
    [~,fName{2}] = system(['find ' appPath filesep appStr filesep ' -type d | grep "/@' name1{1} '"']);
    [~,fName{3}] = system(['find ' appPath filesep appStr filesep ' -type d | grep "/+' name1{1} '"']);
    
    fName = fName(cellfun(@(C)~isempty(C),fName));
    
    if ~isempty(fName)
        fName = strtrim(fName{1});
        % check is there is Contents.m file
        if exist([fName filesep 'Contents.m'],'file')
            helpStr = strtrim(strsplit(fileread(strtrim([fName filesep 'Contents.m'])),'\n'));
            % remove '%' signs
            helpStr(cellfun(@(C)numel(C),helpStr,'uniformoutput',true)==0) = [];
            helpStr = cellfun(@(C)C(2:end),helpStr,'uniformoutput',false);
            helpStr = sprintf('%s\n',helpStr{:});
        end
    end
end


end