function hLine = line(varargin)
% draws a 3D line using patch
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

% number of line segments
nSegment = size(rStart,2);

% create the vertices
V = reshape(permute(cat(3,rStart,rEnd),[1 3 2]),3,[])';

L = (1:nSegment)';
F = [2*L-1 2*L];

% black color
C = repmat([0 0 0],[size(V,1) 1]);

if isempty(hLine)
    % create new patch
    hLine = patch(hAxis,'Vertices',V,'Faces',F,'FaceLighting','flat',...
        'EdgeColor','flat','FaceColor','none','Tag','line','FaceVertexCData',C);
else
    % add to existing patch
    V0 = get(hLine,'Vertices');
    F0 = get(hLine,'Faces');
    C0 = get(hLine,'FaceVertexCData');
    % number of existing faces
    nV0 = size(V0,1);
    set(hLine,'Vertices',[V0;V],'Faces',[F0;F+nV0],'FaceVertexCData',[C0;C]);
end

end