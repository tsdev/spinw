function hText = text(varargin)
% draws a text at a point in 3D
%
% hText = SWPLOT.TEXT(r, string)
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
%
% See also TEXT.
%

if numel(varargin{1}) == 1
    % first input figure handle
    hAxis   = varargin{1};
    r       = varargin{2};
    string  = varargin{3};
else
    hAxis   = gca;
    r       = varargin{1};
    string  = varargin{2};
end

if nargin == 0
    help swplot.text
    return
end

if numel(r) == 3
    r = r(:);
end

nText = size(r,2);

if ~iscell(string)
    string = {string};
end

fontSize = swpref.getpref('fontsize',[]);

hText = gobjects(1,nText);

for ii = 1:nText
    hText(ii) = text(hAxis,r(1,ii),r(2,ii),r(3,ii),string{ii},'FontSize',fontSize,'Color','k',...
    'VerticalAlignment','middle','HorizontalAlignment','center','Tag','text');
end

end