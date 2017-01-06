function handle = text(varargin)
% draws a text at a point in 3D
%
% handle = SWPLOT.TEXT(r, string)
%
% Input:
%
% r         Coordinate of the center of the text.
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

handle = text(hAxis,r(1),r(2),r(3),string,'FontSize',fontSize,'Color','k',...
    'VerticalAlignment','middle','HorizontalAlignment','center','Tag','text');

end