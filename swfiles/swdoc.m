function varargout = swdoc(funName0,varargin)
% opens the SpinW documentation
%
% ### Syntax
%
% `swdoc(funName)`
%
% `link = swdoc(funName)`
%
% ### Description
%
% `swdoc(funName)` shows the documentation on the given function name,
% method, property name. Package functions can be referenced as
% `'packagename.functionname'`, e.g. `'swpref.getpref'`, class methods and
% properties can be referenced the same way, e.g. `'spinw.genmagstr'`. Also
% class, package or folder (that contains `.m` files) names can be
% referenced as well, e.g. `'swsym'`.
%
% `link = swdoc(funName)` returns the web link pointing to the right
% documentation.
%

pref = swpref;

if nargin == 0
    funName0 = '';
end

if strcmp(funName0,'server') && ~isempty(varargin)
    % start/stop Jekyll server
    % stop server if there were any running
    system('pkill jekyll');
    
    switch varargin{1}
        case 'start'
            % start server only works on macOS
            path0 = pwd;
            % TODO change location of documentation files
            if numel(varargin)>1
                cd(varargin{2})
            else
                cd('~/spinwdoc_git')
            end
            system('source ~/.bash_profile; bundle exec jekyll serve --watch&','-echo');
            cd(path0)
            % use the local documentation
            swpref.setpref('docurl','http://localhost:4002')
        case 'stop'
            % already stopped
            % use the online documentation
            swpref.setpref('docurl','https://tsdev.github.io/spinwdoc')
        otherwise
            error('swdoc:WrongParameter','The second input parameter has to be ''start''/''stop'' to start/stop the Jekyll server!')
    end
    varargout = {};
    return
end

funName = funName0;

% list of supported classes where method/property names can be defined
classList = {'spinw'};
idx       = 1;

while ~any(funName=='.') && idx <= numel(classList)
    % check if funName is a SpinW method
    mList = methods(classList{idx});
    pList = properties(classList{idx});
    % find method
    if ~isempty(find(ismember([mList; pList],funName),1)) && ~strcmp(funName,classList{idx})
        funName = [classList{idx} '.' funName]; %#ok<AGROW>
    end
    idx = idx+1;
end

funName(funName=='@') = [];
funName(funName=='.') = '_';

% open documentation using stored url + function name
<<<<<<< Updated upstream
docUrl = pref.docurl;
=======
docUrl = swpref.getpref('docurl',[]);
>>>>>>> Stashed changes
if ~isempty(docUrl) && docUrl(end)~='/'
    docUrl = [docUrl '/'];
end

link = [docUrl funName];

% test if link exists
try
    webread(link);
    if nargout > 0
        varargout{1} = link;
    else
        web(link);
    end
catch
    if isempty(funName0)
        fprintf('Documentation is not found on the given address (%s).\n',swpref.getpref('docurl',[]))
    else
        disp([funName0 ' not found.']);
    end
end

end