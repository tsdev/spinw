function hText = text(varargin)
% creates text at a 3D position
% 
% ### Syntax
% 
% `hText = swplot.text(r, string)`
% 
% `hText = swplot.text(r, string, fontSize)`
%
% `hText = swplot.text(handle, ...)`
%
% ### Description
% 
% `hText = swplot.text(r, string)` creates single or multiple text in 3D
% space.
%  
% `hPatch = swplot.text(handle, ...)` adds the generated text object to a
% given axis referenced by `handle`.
%  
% ### Input Arguments
% 
% `handle`
% : Handle of an axis object, default value is [matlab.gca].
% 
% `r`
% : Coordinate of the center of the text for a single text or
%   matrix with dimensions $[3\times n_{obj}]$ for multiple text.
% 
% `string`
% : String that contains the text or cell of strings when multiple
%   text is drawn.
% 
% `fontSize`
% : Font size in pt, default value is stored in `swpref.getpref('fontsize')`.
% 
% ### See Also
% 
% [matlab.text]
%

if nargin == 0
    swhelp swplot.text
    return
end

fontSize = [];
pref = swpref;

if numel(varargin{1}) == 1
    % first input figure handle
    hAxis   = varargin{1};
    r       = varargin{2};
    string  = varargin{3};
    if nargin>3
        fontSize = varargin{4};
    end
else
    hAxis   = gca;
    r       = varargin{1};
    string  = varargin{2};
    if nargin>3
        fontSize = varargin{3};
    end
end

if numel(r) == 3
    r = r(:);
end

nText = size(r,2);

if ~iscell(string)
    string = {string};
end

if isempty(fontSize)
    fontSize = pref.fontsize;
end

hText = gobjects(1,nText);

for ii = 1:nText
    hText(ii) = text(r(1,ii),r(2,ii),r(3,ii),string{ii},'Parent',hAxis,'FontSize',fontSize,'Color','k',...
    'VerticalAlignment','middle','HorizontalAlignment','center','Tag','text');
end

end