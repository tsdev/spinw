function hPatch = line(varargin)
% draws a 3D line using patch
%
% hLine = SWPLOT.LINE(rStart, rEnd, lineStyle)
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
% lineStyle Line style, default is '-' for continuous line.
%
% See also LINE.
%

if nargin == 0
    help swplot.line
    return
end

lineStyle = '-';
lineWidth = 0.5;

if numel(varargin{1}) == 1
    % first input figure/patch handle
    if strcmp(get(varargin{1},'Type'),'axes')
        hAxis  = varargin{1};
        hPatch = [];
    else
        hAxis  = gca;
        hPatch = varargin{1};
    end
    rStart  = varargin{2};
    rEnd    = varargin{3};
    if nargin > 3
        lineStyle = varargin{4};
    end
    if nargin > 4
        lineWidth = varargin{5};
    end
else
    hAxis   = gca;
    hPatch  = [];
    rStart  = varargin{1};
    rEnd    = varargin{2};
    if nargin > 2
        lineStyle = varargin{3};
    end
    if nargin > 3
        lineWidth = varargin{4};
    end

end

if numel(rStart) == 3
    rStart = rStart(:);
    rEnd   = rEnd(:);
else
    if size(rStart,1)~=3 || size(rEnd,1)~=3
        error('line:WrongInput','To plot multiple line segment use [3 nSegment] matrices!');
    end
end

% number of line segments
nObject = size(rStart,2);

% create the vertices
V = reshape(permute(cat(3,rStart,rEnd),[1 3 2]),3,[])';

L = (1:nObject)';
F = [2*L-1 2*L];

% black color
C = repmat([0 0 0],[size(V,1) 1]);

if strcmp(get(hAxis,'Tag'),'swaxis') && strcmp(lineStyle,'-') && lineWidth == 0.5
    % make sure we are on the plot axis of an swobject
    % add object to the existing edge patch, but only for continuous lines
    hFigure = get(hAxis,'Parent');
    hPatch = getappdata(hFigure,'edgepatch');
end

if isempty(hPatch)
    % create new patch
    hPatch = patch('Parent',hAxis,'Vertices',V,'Faces',F,'FaceLighting','flat',...
        'EdgeColor','flat','FaceColor','none','Tag','line',...
        'FaceVertexCData',C,'LineStyle',lineStyle,'LineWidth',lineWidth);
else
    % add to existing patch
    V0 = get(hPatch,'Vertices');
    F0 = get(hPatch,'Faces');
    C0 = get(hPatch,'FaceVertexCData');
    % number of existing faces
    nV0 = size(V0,1);
    set(hPatch,'Vertices',[V0;V],'Faces',[F0;F+nV0],'FaceVertexCData',[C0;C]);
end

if strcmp(get(hAxis,'Tag'),'swaxis')
    % replicate the arrow handle to give the right number of added objects
    hPatch = repmat(hPatch,[1 nObject]);
end

end