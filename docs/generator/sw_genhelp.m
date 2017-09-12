function doctree = sw_genhelp(varargin)
% generates markdown files from function help
%
% SW_GENHELP('option1', value1, ...)
%

inpForm.fname  = {'path' 'sidebar'   };
inpForm.defval = {{}     'sw_sidebar'};
inpForm.size   = {[1 -1] [1 -2]      };

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
    
    % type of folder
    doctree(ii).isPackage = false;
    doctree(ii).isClass   = false;
    
    switch doctree(ii).name(1)
        case '+'
            % package
            doctree(ii).isPackage = true;
            doctree(ii).fun = doctree(ii).name(2:end);
        case '@'
            % class
            doctree(ii).isClass = true;
            doctree(ii).fun = doctree(ii).name(2:end);
        otherwise
            doctree(ii).fun = doctree(ii).name;
    end

    if ~doctree(ii).isClass
        % find all *.m files in the folder
        fList = dir([path0{ii} filesep '*.m']);
        fList = {fList(:).name};
        doctree(ii).content = struct('file',cell(1,numel(fList)));
        [doctree(ii).content(:).file] = fList{:};
    else
        
        propNames   = properties(doctree(ii).fun);
        methodNames = methods(doctree(ii).fun);
        
        doctree(ii).content = struct('fun',cell(1,numel(allNames)));
        [doctree(ii).content(:).fun] = allNames{:};
    end
    
end

for ii = 1:nPath
    % remove old documents
    try
        rmdir([docroot 'pages' filesep doctree(ii).name],'s')
    catch
    end
    % add new empty folder
    mkdir([docroot 'pages' filesep doctree(ii).name]);
end

for ii = 1:nPath
    % load the help text for each file
    for jj = 1:numel(doctree(ii).content)
        if doctree(ii).isClass
            doctree(ii).content(jj).text = strsplit(help([doctree(ii).fun '.' doctree(ii).content(jj).fun]),newline);
        else
            doctree(ii).content(jj).text = strsplit(help([doctree(ii).folder filesep doctree(ii).content(jj).file]),newline);
        end
        % remove common leading spaces
        doctree(ii).content(jj).text = rmspace(doctree(ii).content(jj).text);
    end
    
    % find Contents.m file and put it to the first place
    if doctree(ii).isClass
        isContents = ismember({doctree(ii).content(:).fun},doctree(ii).fun);
    else
        isContents = ismember({doctree(ii).content(:).file},'Contents.m');
    end
    
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
        doctree(ii).content(jj) = helpcontentfun(doctree(ii).content(jj),doctree(ii).fun,doctree(ii).isPackage);
    end
end

for ii = 1:nPath
    doctree(ii).content(1).fun = '';
    doctree(ii).content(1).frontmatter = struct;
    
    for jj = 1:numel(doctree(ii).content)
        content = doctree(ii).content(jj);
        
        if ~doctree(ii).isClass
            content.fun     = content.file(1:end-2);
        end
        
        % title
        if content.isContents
            if doctree(ii).isPackage
                content.frontmatter.title = ['Package ' doctree(ii).name(2:end)];
            elseif doctree(ii).isClass
                content.frontmatter.title = ['Class ' doctree(ii).name(2:end)];
            else
                content.frontmatter.title = ['Functions in ' doctree(ii).name];
            end
        else
            if doctree(ii).isPackage
                content.frontmatter.title = [doctree(ii).name(2:end) '.' content.fun '( )'];
            elseif doctree(ii).isClass
                content.frontmatter.title = [doctree(ii).name '.' content.fun '( )'];
            else
                content.frontmatter.title = [content.fun '( )'];
            end
        end
        % summary
        content.frontmatter.summary = strtrim(content.text{1});
        % keywords
        content.frontmatter.keywords  = 'sample';
        % sidebar
        content.frontmatter.sidebar   = param.sidebar;
        % permalink
        if content.isContents
            content.frontmatter.permalink = [doctree(ii).fun '.html'];
        else
            content.frontmatter.permalink = [doctree(ii).fun '_' content.fun '.html'];
        end
        % folder
        content.frontmatter.folder    = doctree(ii).name;
        % mathjax
        isMath    = ~isempty(regexp(content.text,'\$\$','once'));
        falsetrue = {'false' 'true'};
        content.frontmatter.mathjax   = falsetrue{isMath+1};
        doctree(ii).content(jj)       = content;
    end
end

% YAML java
javaaddpath(YAML.jarfile);
% Load yaml into java obj
snakeyaml = org.yaml.snakeyaml.Yaml;

% generate the .md files
for ii = 1:nPath
    for jj = 1:numel(doctree(ii).content)
        
        content     = doctree(ii).content(jj);
        %frontmatter = YAML.dump(content.frontmatter)
        
        % Convert to matlab object
        frontmatter = char(snakeyaml.dump(YAML.dump_data(content.frontmatter)));
        % save the text as .md file
        fid = fopen([docroot 'pages' filesep content.frontmatter.folder filesep content.fun '.md'],'w');
        % add newline
        helpText = cellfun(@(C)[C newline],content.text,'UniformOutput',false);
        fprintf(fid,['---' newline frontmatter newline '---' newline helpText{:}]);
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
for ii = 1:nPath
    sidebar.entries.folders(2).folderitems(ii).title  = doctree(ii).name;
    sidebar.entries.folders(2).folderitems(ii).url    = ['/' doctree(ii).fun '.html'];
    sidebar.entries.folders(2).folderitems(ii).output = 'web, pdf';
    
    nContent = find(~[doctree(ii).content.isContents]);
    for jj = 1:numel(nContent)
        content = doctree(ii).content(nContent(jj));
        if doctree(ii).isPackage
            sidebar.entries.folders(2).folderitems(ii).subfolders.title  = [doctree(ii).fun ' package reference'];
        elseif doctree(ii).isClass
            sidebar.entries.folders(2).folderitems(ii).subfolders.title  = [doctree(ii).fun ' class reference'];
        else
            sidebar.entries.folders(2).folderitems(ii).subfolders.title  = [doctree(ii).fun ' folder reference'];
        end
        sidebar.entries.folders(2).folderitems(ii).subfolders.output = 'web, pdf';
        sidebar.entries.folders(2).folderitems(ii).subfolders.subfolderitems(jj).title  = content.frontmatter.title(1:end-3);
        sidebar.entries.folders(2).folderitems(ii).subfolders.subfolderitems(jj).url    = ['/' content.frontmatter.permalink];
        sidebar.entries.folders(2).folderitems(ii).subfolders.subfolderitems(jj).output = 'web, pdf';
    end
end

%yamlStr = YAML.dump(sidebar);
yamlStr = char(snakeyaml.dump(YAML.dump_data(sidebar)));

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
        ll2 = [folder '.'];
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
            text{ii} = [ll1 '* [' ll2 line{1} '](/' folder '_' line{1} ') ' line{2}];
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