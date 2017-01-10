function hLine = line(varargin)
% draws a 3D line
%
% hLine = SWPLOT.LINE(rStart, rEnd)
%
% hLine = SWPLOT.LINE(handle,...)
%
% Handle can be the handle of an axes object or a line object. It either
% selects an axis to plot or a patch object (triangulated) to add vertices
% and faces.
%
% Input:
%
% handle    Handle of an axis or patch object. In case of patch object, the
%           constructed faces will be added to the existing object instead
%           of creating a new one.
% rStart    Coordinate(s) of the starting point, either a 3 element vector or
%           a matrix with dimensions [3 nLineSegment] to plot multiple line
%           segments.
% rEnd      Coordinate(s) of the end point, either a 3 element vector or
%           a matrix with dimensions [3 nLineSegment] to plot multiple line
%           segments.
%
% See also LINE.
%

if nargin == 0
    help swplot.line
    return
end

if numel(varargin{1}) == 1
    % first input figure/patch handle
    if strcmp(get(varargin{1},'Type'),'axes')
        hAxis  = varargin{1};
        hLine = [];
    else
        hAxis  = gca;
        hLine = varargin{1};
    end
    rStart  = varargin{2};
    rEnd    = varargin{3};
else
    hAxis   = gca;
    hLine  = [];
    rStart  = varargin{1};
    rEnd    = varargin{2};
end

if numel(rStart) == 3
    rStart = rStart(:);
    rEnd   = rEnd(:);
else
    if size(rStart,1)~=3 || size(rEnd,1)~=3
        error('line:WrongInput','To plot multiple line segment use [3 nSegment] matrices!');
    end
end

% create the segments and remove the last nan
r = reshape(permute(cat(3,rStart,rEnd,nan(size(rStart))),[1 3 2]),3,[]);
r = r(:,1:(end-1));

if isempty(hLine)
    % create new patch
    hLine = line(hAxis,r(1,:),r(2,:),r(3,:),'Color','k','LineStyle','-','Tag','line');
else
    r0 = [get(hLine,'XData');get(hLine,'YData');get(hLine,'ZData')];
    
    if any(isnan(r0(:,end)))
        r = [r0 nan(3,1) r];
    else
        r = [r0 r];
    end
    
    % add line segments
    set(hLine,'XData',r(1,:),'Ydata',r(2,:),'ZData',r(3,:));
end

end