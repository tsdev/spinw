function hText = text(varargin)
% draws a text at a point in 3D
%
% hText = SWPLOT.TEXT(r, string, {fontSize})
%
% hText = SWPLOT.TEXT(handle,...)
%
% Handle of an axes object that selects an axis to plot.
%
% Input:
%
% handle    Handle of an axis object.
% r         Coordinate of the center of the text for a single text or
%           matrix with dimensions [3 nText] for multiple text.
% string    String contains the text or cell of strings to plot multiple
%           text.
% fontSize  Font size in pt, default value is stored in
%           swpref.getpref('fontsize')
%
% See also TEXT.
%

if nargin == 0
    help swplot.text
    return
end

fontSize = [];

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
    fontSize = swpref.getpref('fontsize',[]);
end

hText = gobjects(1,nText);

for ii = 1:nText
    hText(ii) = text(r(1,ii),r(2,ii),r(3,ii),string{ii},'Parent',hAxis,'FontSize',fontSize,'Color','k',...
    'VerticalAlignment','middle','HorizontalAlignment','center','Tag','text');
end

end