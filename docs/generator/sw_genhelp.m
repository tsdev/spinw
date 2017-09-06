function sw_genhelp(varargin)
% generates markdown files from function help
%
% SW_GENHELP('option1', value1, ...)
%

inpForm.fname  = {'path'};
inpForm.defval = {{}    };
inpForm.size   = {[1 -1]};

param = sw_readparam(inpForm, varargin{:});

if ~iscell(param.path)
    path0 = {param.path};
else
    path0 = param.path;
end

nPath = numel(path0);

docroot  = [sw_rootdir 'docs' filesep];
doctree  = struct('name',cell(1,nPath),'folder',[],'content',[]);

% loop over all path to generate help files
for ii = 1:nPath
    % name of the parent folder
    [~,pp1,pp2] = fileparts(path0{ii});
    doctree(ii).name   = [pp1 pp2];
    doctree(ii).folder = path0{ii};
    
    % find all *.m files in the folder
    fList = dir([path0{ii} filesep '*.m']);
    fList = {fList(:).name};
    doctree(ii).content = struct('file',cell(1,numel(fList)));
    [doctree(ii).content(:).file] = fList{:};
    
    % type of folder
    doctree(ii).isPackage = false;
    doctree(ii).isClass   = false;
    
    switch doctree(ii).name
        case '+'
            % package
            doctree(ii).isPackage = true;
        case '@'
            % class
            doctree(ii).isClass = true;
    end
    
end

for ii = 1:nPath
    % remove old documents
    try
        rmdir([docroot 'pages' filesep doctree(ii).folder],'s')
    catch
    end
    % add new empty folder
    mkdir([docroot 'pages' filesep doctree(ii).folder]);
end

for ii = 1:nPath
    % load the help text for each file
    for jj = 1:numel(doctree(ii).content)
        doctree(ii).content(jj).text = strsplit(help([doctree(ii).folder filesep doctree(ii).content(jj).file]),newline);
        % remove common leading spaces
        doctree(ii).content(jj).text = rmspace(doctree(ii).content(jj).text);
    end
    
    % find Contents.m file and put it to the first place
    isContents = ismember({doctree(ii).content(:).file},'Contents.m');
    
    isContents = num2cell(isContents);
    [doctree(ii).content(:).isContents] = isContents{:};
    
    % put Contents.m file to the first place
    [~,idx] = sort(cell2mat(isContents),'descend');
    doctree(ii).content = doctree(ii).content(idx);
end

% convert Contents.m files
for ii = 1:nPath
    cIdx = find([doctree(ii).content(:).isContents]);
    for jj = 1:numel(cIdx)
        doctree(ii).content(jj) = helpcontentfun(doctree(ii).content(jj),doctree(ii).name,doctree(ii).isPackage);
    end
end

for ii = 1:nPath
    for jj = 1:numel(doctree(ii).content)
        
        % summary
        doctree(ii).content(jj).summary = strtrim(doctree(ii).content(jj).text{1});
        doctree(ii).content(jj).fun     = doctree(ii).content(jj).file(1:end-2);
        
        % title
        if doctree(ii).content(jj).isContents
            if doctree(ii).isPackage
                doctree(ii).content(jj).title = ['Package ' doctree(ii).name(2:end)];
            elseif doctree(ii).isClass
                doctree(ii).content(jj).title = ['Class ' doctree(ii).name(2:end)];
            else
                doctree(ii).content(jj).title = ['Functions in ' doctree(ii).name];
            end
        else
            if doctree(ii).isPackage
                doctree(ii).content(jj).title = [doctree(ii).name(2:end) '.' doctree(ii).content(jj).fun '( )'];
            elseif doctree(ii).isClass
                doctree(ii).content(jj).title = [doctree(ii).name '.' doctree(ii).content(jj).fun '( )'];
            else
                doctree(ii).content(jj).title = [doctree(ii).content(jj).fun '( )'];
            end
        end
        % permalink
        
        
        frontmatter = YAML.dump(struct(...
            'title',    titleStr0,...
            'keywords', 'sample',...
            'summary',  summary,...
            'sidebar',  'sw_sidebar',...
            'permalink',[folder{ii} '_' fN '.html'],...
            'folder',   folder{ii},...
            'mathjax',  'true'));
        header = ['---' newline frontmatter newline '---' newline];
        
        % save the text as .md file
        fid = fopen([docroot 'pages' filesep folder{ii} filesep folder{ii} '_' fN '.md'],'w');
        % add newline
        helpText = cellfun(@(C)[C newline],helpText,'UniformOutput',false);
        fprintf(fid,[header helpText{:}]);
        fclose(fid);
    end
end

% generate sidebar YAML file
swver = sw_version;
if isempty(swver.Version)
    verStr = ['R' swver.Revision];
else
    verStr = swver.Version;
end

sidebar = struct;
sidebar.entries.title   = 'Sidebar';
sidebar.entries.product = 'SpinW';
sidebar.entries.version = verStr;
% documentation
sidebar.entries.folders(1).title  = 'Documentation';
sidebar.entries.folders(1).output = 'web, pdf';

sidebar.entries.folders(2).title  = 'Function reference';
sidebar.entries.folders(2).output = 'web, pdf';
%sidebar.entries.folders(1).type   = 'frontmatter';
for ii = 1:numel(folder)
    sidebar.entries.folders(2).folderitems(ii).title  = folder{ii};
    sidebar.entries.folders(2).folderitems(ii).url    = ['/' folder{ii} '.html'];
    sidebar.entries.folders(2).folderitems(ii).output = 'web, pdf';
    
    for jj = 1:numel(funName{ii})
        sidebar.entries.folders(2).folderitems(ii).subfolders.title  = ['Functions in ' folder{ii}];
        sidebar.entries.folders(2).folderitems(ii).subfolders.output = 'web, pdf';
        sidebar.entries.folders(2).folderitems(ii).subfolders.subfolderitems(jj).title  = titleStr{ii}{jj};
        sidebar.entries.folders(2).folderitems(ii).subfolders.subfolderitems(jj).url    = ['/' folder{ii} '_' funName{ii}{jj} '.html'];
        sidebar.entries.folders(2).folderitems(ii).subfolders.subfolderitems(jj).output = 'web, pdf';
    end
end

yamlStr = YAML.dump(sidebar);

% add extra '-' that is required by Jekyll
yamlStr1 = strsplit(yamlStr,newline);
nSpace   = zeros(1,numel(yamlStr1));
for ii = 1:numel(yamlStr1)
    temp = find(diff(yamlStr1{ii}~=' '),1,'first');
    if isempty(temp)
        nSpace(ii) = 0;
    else
        nSpace(ii) = temp;
    end
end

iLine  = find(diff(nSpace))+1;
nSpace = nSpace(iLine);

iLine(nSpace==0)  = [];
nSpace(nSpace==0) = [];

for ii = 1:numel(iLine)
    if yamlStr1{iLine(ii)}(nSpace(ii)+1) ~= '-' && yamlStr1{iLine(ii)-1}(nSpace(ii)-1) ~= '-'
        yamlStr1{iLine(ii)}(nSpace(ii)-1) = '-';
    end
end
yamlStr1 = cellfun(@(C)[C newline],yamlStr1,'uniformoutput',false);
yamlStr1 = [yamlStr1{:}];

fid = fopen([docroot filesep '_data' filesep 'sidebars' filesep 'sw_sidebar.yml'],'w');
fprintf(fid,yamlStr1);
fclose(fid);

end

function str = rmspace(str)
% remove leading spaces

% remove empty lines
str(cellfun(@(C)isempty(C),str)) = [];
% lines that begin with space
sIdx = cellfun(@(C)C(1)==' ' && any(C~=' '),str);
% minimum space
nSpace = min(cellfun(@(C)find(diff(C==' '),1,'first'),str(sIdx)));
% remove leading spaces
str(sIdx) = cellfun(@(C)C(nSpace+1:end),str(sIdx),'UniformOutput',false);

end

function doccontent = helpcontentfun(doccontent,folder,isPackage)
% convert Contents.m text

text = doccontent.text;

% find line "Files"
lIdx = find(ismember(doccontent.text,'Files'));
if ~isempty(lIdx)
    text{lIdx} = '### Files';
    
    firstLine = true;
    
    if isPackage
        ll2 = [folder(2:end) '.'];
    else
        ll2 = [];
    end
    
    for ii = lIdx+1:numel(text)
        line = strtrim(strsplit(text{ii},'-'));
        if numel(line)==2 && all(text{ii}~=':')
            % add link
            if firstLine
                ll1 = newline;
            else
                ll1 = [];
            end
            text{ii} = [ll1 '* [' ll2 line{1} '()](/' folder '_' line{1} ') ' line{2}];
            firstLine = false;
        elseif any(text{ii}==':')
            % sub header line
            text{ii} = ['#### ' text{ii}];
            firstLine = true;
        end
    end
    
    doccontent.text = text;
end

end