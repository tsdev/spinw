function zoom(mode, hFigure)
% zooming figure
%
% SWPLOT.ZOOM(mode)
%
% Input:
%
% mode      Either a number determining the zoom value, or 'auto' that
%           zooms to see every object on the figure.
%
% See also SWPLOT.FIGURE.
%

if nargin == 0
    mode = 'auto';
end

if nargin < 2
    % find active figure
    hFigure = swplot.activefigure;
end

hAxis = getappdata(hFigure,'axis');

if isnumeric(mode)    
    set(hAxis,'CameraViewAngle',mode);
elseif strcmpi(mode,'auto')
    set(hAxis,'CameraViewAngleMode','auto');
    set(hAxis,'CameraViewAngle',2*get(hAxis,'CameraViewAngle'));
else
    error('zoom:WrongInput','Wrong zoom mode!');
end

end