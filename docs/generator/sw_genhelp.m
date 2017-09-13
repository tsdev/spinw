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
    doctree(ii).fullname = [pp1 pp2];
    doctree(ii).folder   = path0{ii};
    
    % type of folder
    doctree(ii).isPackage = false;
    doctree(ii).isClass   = false;
    
    switch doctree(ii).fullname(1)
        case '+'
            % package
            doctree(ii).isPackage = true;
            doctree(ii).name      = doctree(ii).fullname(2:end);
        case '@'
            % class
            doctree(ii).isClass = true;
            doctree(ii).name = doctree(ii).fullname(2:end);
        otherwise
            doctree(ii).name = doctree(ii).fullname;
    end

    if doctree(ii).isClass
        % class
        name        = doctree(ii).name;
        
        propNames   = properties(name);
        methodNames = methods(name);
        
        funList = [name; cellfun(@(C)[name '.' C],[propNames;methodNames],'UniformOutput',false)];
        doctree(ii).content = struct('fun',cell(1,numel(funList)),'isProp',[],'file',[]);
        [doctree(ii).content(:).fun] = funList{:};
        isProp = num2cell([false true(1,numel(propNames)) false(1,numel(methodNames))]);
        [doctree(ii).content(:).isProp] = isProp{:};
    elseif doctree(ii).isPackage
        % package
        name        = doctree(ii).name;
        
        % find all *.m files in the folder
        fList = dir([path0{ii} filesep '*.m']);
        fList = {fList(:).name};
        % remove .m
        nList = cellfun(@(C)[name '.' C(1:end-2)],fList,'UniformOutput',false);
        doctree(ii).content = struct('file',cell(1,numel(fList)),'fun',[]);
        [doctree(ii).content(:).file] = fList{:};
        [doctree(ii).content(:).fun]  = nList{:};
    else
        % find all *.m files in the folder
        fList = dir([path0{ii} filesep '*.m']);
        fList = {fList(:).name};
        % remove .m
        nList = cellfun(@(C)C(1:end-2),fList,'UniformOutput',false);
        doctree(ii).content = struct('file',cell(1,numel(fList)),'fun',[]);
        [doctree(ii).content(:).file] = fList{:};
        [doctree(ii).content(:).fun]  = nList{:};
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
        if isempty(doctree(ii).content(ii).file)
            doctree(ii).content(jj).text = strsplit(help(doctree(ii).content(jj).fun),newline);
        else
            doctree(ii).content(jj).text = strsplit(help([doctree(ii).folder filesep doctree(ii).content(jj).file]),newline);
        end
        % remove common leading spaces
        doctree(ii).content(jj).text = rmspace(doctree(ii).content(jj).text);
    end
    
    % find Contents.m file and put it to the first place
    if doctree(ii).isClass
        isContents = ismember({doctree(ii).content(:).fun},doctree(ii).name);
    else
        isContents = ismember({doctree(ii).content(:).file},'Contents.m');
    end
    
    isContents = num2cell(isContents);
    [doctree(ii).content(:).isContents] = isContents{:};
    
    % put Contents.m file to the first place
    [~,idx] = sort(cell2mat(isContents),'descend');
    doctree(ii).content = doctree(ii).content(idx);
end

% convert Contents.m files and mine out titles
for ii = 1:nPath
    cIdx = find([doctree(ii).content(:).isContents]);
    for jj = 1:numel(cIdx)
        [doctree(ii).content(cIdx(jj)), tStruct] = helpcontentfun(doctree(ii).content(cIdx(jj)),doctree(ii).name,doctree(ii).isPackage || doctree(ii).isClass);
    end
    % assign title from tList
    fList  = {doctree(ii).content(:).fun};
        
    if all(cellfun(@(C)isempty(C),{tStruct.title}))
        % no titles are defined
        if doctree(ii).isClass
            title0 = 'Methods';
        else
            title0 = 'Files';
        end
        title0 = repmat({title0},1,numel(tStruct));
        [tStruct(:).title] = title0{:};
    end
    
    doctree(ii).utitle = ['Properties' unique({tStruct.title},'stable')];
    % add properties
    if doctree(ii).isClass && any(~[doctree(ii).content.isProp])
        pList = {doctree(ii).content([doctree(ii).content.isProp]).fun};
        [tStruct(end+(1:numel(pList))).fun] = pList{:};
        title0 = repmat({'Properties'},1,numel(pList));
        [tStruct(end-(numel(pList):-1:1)+1).title] = title0{:};
    end
    
    tfList = {tStruct(:).fun};
    
    for jj = 1:numel(fList)
        idx = find(ismember(tfList,fList{jj}));
        if numel(idx) == 1
            doctree(ii).content(jj).title = tStruct(idx).title;
        else
            doctree(ii).content(jj).title = 'Miscellaneous';
            if ~doctree(ii).content(jj).isContents && ~strcmp(doctree(ii).content(jj).fun,doctree(ii).name)
                warning([doctree(ii).content(jj).fun ' has no corresponding title defined, please append Contents.m or class description!'])
            end
        end
    end
    
end

for ii = 1:nPath
    doctree(ii).content(1).frontmatter = struct;
    
    for jj = 1:numel(doctree(ii).content)
        content = doctree(ii).content(jj);
        
        % title
        if content.isContents
            if doctree(ii).isPackage
                content.frontmatter.title = ['Package ' doctree(ii).name];
            elseif doctree(ii).isClass
                content.frontmatter.title = ['Class ' doctree(ii).name];
            else
                content.frontmatter.title = ['Functions in ' doctree(ii).name];
            end
        else
            content.frontmatter.title = [content.fun '( )'];
        end
        % summary
        content.frontmatter.summary = strtrim(content.text{1});
        % keywords
        content.frontmatter.keywords  = 'sample';
        % sidebar
        content.frontmatter.sidebar   = param.sidebar;
        % permalink
        if content.isContents
            content.frontmatter.permalink = [doctree(ii).name '.html'];
        else
            content.frontmatter.permalink = [strrep(content.fun,'.','_'),'.html'];
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

for ii = 1:nPath
    sidebar.entries.folders(ii+1).title  = doctree(ii).name;
    sidebar.entries.folders(ii+1).output = 'web, pdf';


    sidebar.entries.folders(ii+1).folderitems(1).title  = 'Description';
    sidebar.entries.folders(ii+1).folderitems(1).url    = ['/' doctree(ii).name '.html'];
    sidebar.entries.folders(ii+1).folderitems(1).output = 'web, pdf';
    
    % find unique titles
    tCell   = {doctree(ii).content.title};
    %tUnique = unique(tCell,'stable');
    tUnique = doctree(ii).utitle;
    mIdx = find(ismember(tUnique,'Miscellaneous'));
    if ~isempty(mIdx)
        tUnique = tUnique([1:(mIdx-1) (mIdx+1):end mIdx]);
    end
    
    for jj = 1:numel(tUnique)
        idx = find(ismember(tCell,tUnique{jj}) & ~[doctree(ii).content.isContents]);
        if ~isempty(idx)
            % title
            sidebar.entries.folders(ii+1).folderitems(1).subfolders(jj).title = tUnique{jj};
            sidebar.entries.folders(ii+1).folderitems(1).subfolders(jj).output = 'web, pdf';
            
            for kk = 1:numel(idx)
                content = doctree(ii).content(idx(kk));
                sidebar.entries.folders(ii+1).folderitems(1).subfolders(jj).subfolderitems(kk).title  = content.frontmatter.title(1:end-3);
                sidebar.entries.folders(ii+1).folderitems(1).subfolders(jj).subfolderitems(kk).url    = ['/' content.frontmatter.permalink];
                sidebar.entries.folders(ii+1).folderitems(1).subfolders(jj).subfolderitems(kk).output = 'web, pdf';
            end
        end
    end
    
end

%yamlStr = YAML.dump(sidebar);
yamlStr  = char(snakeyaml.dump(YAML.dump_data(sidebar)));
yamlStr1 = strsplit(yamlStr,newline);

% add extra newline where {} unit is inline with some previous text
bLine = regexp(yamlStr1,': {');
bLine2 = find(cellfun(@(C)~isempty(C),bLine));
bLine  = bLine(bLine2);
for ii = 1:numel(bLine2)
    newLine = yamlStr1{bLine2(ii)}(bLine{ii}+2:end);
    yamlStr1{bLine2(ii)} = yamlStr1{bLine2(ii)}(1:bLine{ii});
    yamlStr1 = yamlStr1([1:bLine2(ii) bLine2(ii) (bLine2(ii)+1):end]);
    nSpace = sum(yamlStr1{bLine2(ii)}==' ');
    yamlStr1{bLine2(ii)+1} = [repmat(' ',1,nSpace) '- ' newLine];
    bLine2(ii+1:end) = bLine2(ii+1:end)+1;
end

% add extra '-' that is required by Jekyll
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

function [doccontent, tList] = helpcontentfun(doccontent,name,isDot)
% convert Contents.m text

text = doccontent.text;

oldheader = {'Files:' [name ' methods:']};
newheader = {'### Files' '### Methods'};
% find line "Files" / "Methods"
[lIdx, tIdx] = ismember(doccontent.text,oldheader);
lIdx = find(lIdx);
if numel(lIdx) == 1
    
    tList = struct('fun',{},'title',{});
    tIdx = tIdx(lIdx);
    text{lIdx} = newheader{tIdx};
    
    firstLine = true;
    
    if isDot
        ll2 = [name '.'];
    else
        ll2 = [];
    end
    
    title = '';
    
    for ii = lIdx+1:numel(text)
        mIdx = find(text{ii}=='-',1,'first');
        if ~isempty(mIdx) && all(text{ii}~=':')
            line = strtrim({text{ii}(1:mIdx-1) text{ii}(mIdx+1:end)});
            % add link
            if firstLine
                ll1 = newline;
            else
                ll1 = [];
            end
            text{ii} = [ll1 '* [' ll2 line{1} '](/' name '_' line{1} ') ' line{2:end}];
            firstLine = false;
            tList(end+1).fun = [ll2 line{1}]; %#ok<AGROW>
            tList(end).title = title;
        elseif any(text{ii}==':')
            % sub header line
            title = strtrim(text{ii});
            title = title(1:end-1);
            text{ii} = ['#### ' text{ii}];
            firstLine = true;
        end
        
    end
    
    doccontent.text = text;
else
    error('sw_genhelp:MissingContent','Contents.m file is missing for %s!',name);
end

end