function helpStr = sw_help(name, hRoot, showFigure)
% provide help string in a compiled application
%
% SW_HELP(name, path, show)
%
% Input:
%
% name      Name of function, class method, package of folder.
% path      Path, where the function is searching for matching files,
%           default value is the current folder. In deployed code the
%           default is the location of the compiled excutable.
% show      If true, the help text is shown in a new window, default is
%           false.
%

helpStr = '';

if nargin<3
    showFigure = false;
end

if ~isdeployed && nargin < 2
    hRoot = pwd;
end

if nargin < 2
    % get the location of the m-files
    [appPath, appName] = sw_apppath;
    appStr = [ appName filesep 'Source'];
    hRoot  = [appPath filesep appStr];
end

% name without .m extension
name = strsplit(name,'.');

ismethod = false;

if numel(name) == 2 && numel(name{2}) == 1 && name{2}=='m'
    % m-file was given
    sname{1}   = ['/' name{1} '\.m'];
    funName{1} = [name{1} '()'];
elseif numel(name) == 2
    % class method or packaged file
    sname{1} = ['/@' name{1} '/' name{2} '\.m'];
    sname{2} = ['/+' name{1} '/' name{2} '\.m'];
    funName{1} = ['@' name{1} '.' name{2} '()'];
    funName{2} = ['+' name{1} '.' name{2} '()'];
    ismethod = true;
elseif numel(name) == 1
    % m-file file or package or folder name
    sname{1}   = ['/'  name{1} '\.m'];
    funName{1} = [name{1} '()'];
    sname{2}   = ['/+' name{1} '/' 'Contents\.m'];
    funName{2} = ['+' name{1}];
    sname{3}   = ['/'  name{1} '/' 'Contents\.m'];
    funName{3} = [name{1}];
else
    return
end

fName = cell(1,numel(sname));
for ii = 1:numel(sname)
    % find the right file
    [~,fName{ii}] = system(['find ' hRoot ' -type f | grep "' sname{ii} '"']);
end

% keep only non-empty strings
fIdx = cellfun(@(C)~isempty(C),fName);
funName = funName(fIdx);
fName   = fName(fIdx);

if ~isempty(fName)
    % check if mutiple files were found, in that case take the first
    fName = strsplit(fName{1},'\n');
    hFile = strsplit(fileread(fName{1}),'\n');
    
elseif ismethod
    % search for function names inside class definition file
    sname = ['/@' name{1} '/' name{1} '\.m'];
    [~,fName] = system(['find ' hRoot ' -type f | grep "' sname '"']);
    
    if isempty(fName)
        return
    end
    fName   = strsplit(fName,'\n');
    hFile = strsplit(fileread(fName{1}),'\n');
    % find function definition
    lIdx = find(cellfun(@(C)~isempty(C),regexp(hFile,['\s*function.+?=\s*' name{2}])));
    if isempty(lIdx)
        return
    else
        hFile = hFile(lIdx+1:end);
    end
else
    return
end

% extract help text from .m files (first continuous block of comments)
% remove whitespaces
hFile = strtrim(hFile);

% remove the line that defines a function if it is before the first
% commented line
fIdx = find(cellfun(@(C)~isempty(C),regexp(hFile,'\s*(function|classdef)\s+')),1,'first');
cIdx = find(cellfun(@(C)~isempty(C),regexp(hFile,'^%')),1,'first');
if ~isempty(fIdx) && cIdx>fIdx
    hFile(fIdx) = [];
end

% find the first command line and remove all following text
cIdx = find(cellfun(@(C)~isempty(C),regexp(hFile,'^[a-zA-Z]+')),1,'first');
if ~isempty(cIdx)
    hFile = hFile(1:cIdx-1);
end

% remove empty lines
hFile(cellfun(@(C)numel(C),hFile)==0) = [];
% find commented lines
isComment = cellfun(@(C)C(1)=='%',hFile);
% first commented line
lIdx1 = find(isComment,1,'first');
% last commented line
lIdx2 = find(diff(isComment)==-1,1,'first');
if isempty(lIdx1)
    return
end
if isempty(lIdx2)
    lIdx2 = numel(hFile);
end
% keep the commented lines only
hFile = cellfun(@(C)C(2:end),hFile(lIdx1:lIdx2),'uniformoutput',false);
% create string from cell
helpStr = sprintf('%s\n',hFile{1:end-1});
helpStr = [helpStr hFile{end}];

if showFigure
    % show help text in new window
    hFigure = figure('MenuBar','none','DockControls','off','name',funName{1});
    panel1 = uipanel('Parent',hFigure);
    panel2 = uipanel('Parent',panel1);
    set(panel1,'Position',[0 0 1 1]);
    set(panel2,'Position',[0 -1 1 5],'backgroundcolor','w');
    hAxis = axes('Parent',panel2,'color','w');
    set(hAxis,'Position',[0 0 1 1],'clipping','off');
    axis('off');
    text(0,1,helpStr,'Units','normalized','VerticalAlignment','top','interpreter','none');
    uicontrol('Style','Slider','Parent',hFigure,...
        'Units','normalized','Position',[0.95 0 0.05 1],...
        'Value',1,'Callback',{@slider_callback,panel2});
end

end

function slider_callback(src,~,arg1)
% slider callback

val = get(src,'Value');
set(arg1,'Position',[0 -val 1 2])
end