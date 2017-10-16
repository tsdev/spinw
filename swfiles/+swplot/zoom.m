function zoom(mode, hFigure)
% zooms to objects
% 
% ### Syntax
% 
% `swplot.zoom(mode)`
% 
% `swplot.zoom(mode, hFigure)`
%
% ### Description
% 
% `swplot.zoom(mode)` controls the zoom (angle of view of the virtual
% camera) level on the active [swplot] figure.
% 
% `swplot.zoom(mode, hFigure)` controls the zoom on the swplot figure
% referenced by the `hFigure` handle.
%
% ### Input Arguments
% 
% `mode`
% : Either a number determining the relative zoom value, or `'auto'`
%   that zooms to fit every object into the figure.
% 
% `hFigure`
% : Handle of the swplot figure, default value is the active swplot figure.
% 
% ### See Also
% 
% [swplot.figure]
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
    
    % find outer positions of objects
    obj = getappdata(hFigure,'objects');
    % position in lu
    if ~isempty(obj)
        pos = cat(3,obj(:).position);
        size = diff([min(pos(:,1,:),[],3) max(pos(:,1,:),[],3)],[],2);
        size = max(swplot.base(hFigure)*size);
        
        if size>0
            % view angle
            angle = atand(size/norm(get(hAxis,'CameraPosition')));
            
            set(hAxis,'CameraViewAngleMode','auto');
            %set(hAxis,'CameraViewAngle',6*get(hAxis,'CameraViewAngle'));
            set(hAxis,'CameraViewAngle',1.5*angle);
        end
    end
else
    error('zoom:WrongInput','Wrong zoom mode!');
end

end