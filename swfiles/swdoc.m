function varargout = swdoc(funName0)
% opens the SpinW help
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

if nargin == 0
    funName0 = '';
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

if isempty(swpref.getpref('docport',true))
    % read documentation from the web
    link = [swpref.getpref('doclink',true) funName];
else
    % read documentation from localhost with the given port number
    link = ['http://localhost:' num2str(swpref.getpref('docport',true)) '/' funName];
end

% test if link exists
try
    webread(link);
    if nargout > 0
        varargout{1} = link;
    else
        %web(link,'-notoolbar');
        web(link);
    end
catch
    disp([funName0 ' not found.']);
end

end