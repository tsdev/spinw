function hPatch = ellipsoid(varargin)
% draw ellipsoid
%
% hPatch = SWPLOT.ELLIPSOID(R0,T,mesh)
%
% The function can draw multiple ellipsoids with a single patch command.
% Significant speedup can be achieved by a single patch compared to ellipse
% per patch.
%
% hPatch = SWPLOT.ELLIPSOID(handle,...
%
% Handle can be the handle of an axes object or a patch object. It either
% selects an axis to plot or a patch object (triangulated) to add vertices
% and faces.
%
% Input:
%
% handle    Handle of an axis or patch object. In case of patch object, the
%           constructed faces will be added to the existing object instead
%           of creating a new one.
% R0        Center of the ellipsoid stored in a matrix with dimensions of
%           [3 nEllipse].
% T         Transformation matrix that transforms a unit sphere to the
%           ellipse via: R' = T(:,:,i)*R
%           Dimensions are [3 3 nEllipse].
% mesh      Mesh of the ellipse surface, a triangulation class object or an
%           integer that used to generate an icosahedron mesh with #mesh
%           number of additional triangulation.
%
% See also TRIANGULATION, SWPLOT.ICOMESH.
%

if nargin == 0
    help swplot.ellipsoid
    return
end

if numel(varargin{1}) == 1
    % first input figure/patch handle
    if strcmp(get(varargin{1},'Type'),'axes')
        hAxis  = varargin{1};
        hPatch = [];
    else
        hAxis  = gca;
        hPatch = varargin{1};
    end
    
    R0  = varargin{2};
    T    = varargin{3};
    if nargin > 3
        mesh = varargin{4};
    else
        mesh = [];
    end
else
    hAxis  = gca;
    hPatch = [];
    R0  = varargin{1};
    T    = varargin{2};
    if nargin > 2
        mesh = varargin{3};
    else
        mesh = [];
    end
end

if isempty(mesh)
    mesh = 1;
end

if isnumeric(mesh)
    % limit the largest mesh to plot to avoid slowing down Matlab too much
    mesh = min(mesh,swpref.getpref('maxmesh',[]));
    % generate mesh
    mesh = swplot.icomesh(mesh);
end

if numel(R0) == 3
    R0 = R0(:);
end

% plot multiple ellipse
nEllipse = size(R0,2);

sT = [size(T) 1];
if size(R0,1)~=3 || any(sT(1:3)~=[3 3 nEllipse])
    error('ellipsoid:WrongInput','Matrices have incompatible dimensions!')
end

V0 = mesh.Points;
F0 = mesh.ConnectivityList;
NV = size(V0,1);
NF = size(F0,1);

V = zeros(NV*nEllipse,3);
F = zeros(NF*nEllipse,3);

for ii = 1:nEllipse
    % vertices
    V((1:NV)+(ii-1)*NV,:) = bsxfun(@plus,T(:,:,ii)*V0',R0(:,ii))';
    % faces
    F((1:NF)+(ii-1)*NF,:) = F0+(ii-1)*NV;
end

% red color
C = repmat([1 0 0],[size(F,1) 1]);

if isempty(hPatch)
    % create patch
    hPatch = patch(hAxis,'Vertices',V,'Faces',F,'FaceLighting','flat',...
        'EdgeColor','none','FaceColor','flat','Tag','ellipsoid','FaceVertexCData',C);
else
    V0 = get(hPatch,'Vertices');
    F0 = get(hPatch,'Faces');
    C0 = get(hPatch,'FaceVertexCData');
    % number of existing faces
    nV0 = size(V0,1);
    set(hPatch,'Vertices',[V0;V],'Faces',[F0;F+nV0],'FaceVertexCData',[C0;C]);
end

end