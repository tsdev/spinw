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
% r         Coordinate of the center of the text for a single text.
% string    String contains the text.
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

fontSize = swpref.getpref('fontsize',[]);

hText = text(hAxis,r(1),r(2),r(3),string,'FontSize',fontSize,'Color','k',...
    'VerticalAlignment','middle','HorizontalAlignment','center','Tag','text');

end