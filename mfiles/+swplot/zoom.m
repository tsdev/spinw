function zoom(mode, hFigure)
% zooming figure
%
% SWPLOT.ZOOM(mode, {hFigure})
%
% Input:
%
% mode      Either a number determining the relative zoom value, or 'auto'
%           that zooms to see every object on the figure.
% hFigure   Handle of the swplot figure window, optional.
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
    set(hAxis,'CameraViewAngle',get(hAxis,'CameraViewAngle')/mode);
elseif strcmpi(mode,'auto')
    set(hAxis,'CameraViewAngleMode','auto');
    set(hAxis,'CameraViewAngle',1.5*get(hAxis,'CameraViewAngle'));
else
    error('zoom:WrongInput','Wrong zoom mode!');
end

end