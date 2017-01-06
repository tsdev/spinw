function handle = line(varargin)
% draws a 3D line
%
% handle = SWPLOT.LINE(rStart, rEnd)
%
% handle = SWPLOT.LINE(hAxis,...)
%
% Plots to specific axis.
%
% Input:
%
% hAxis     Axis handle, default is gca.
% rStart    Coordinate of the starting point.
% rEnd      Coordinate of the end point.
%
% See also LINE.
%

if nargin == 0
    help swplot.line
    return
end

if numel(varargin{1}) == 1
    % first input figure handle
    hAxis   = varargin{1};
    rStart  = varargin{2};
    rEnd    = varargin{3};
else
    hAxis   = gca;
    rStart  = varargin{1};
    rEnd    = varargin{2};
end

r = [rStart(:) rEnd(:)];
handle = line(hAxis,r(1,:),r(2,:),r(3,:),'Color','k','LineStyle','-','Tag','line');

end