function handle = ellipsoid(varargin)
% draw ellipsoid
%
% handle = SWPLOT.ELLIPSOID(R0,T,mesh)
%
% The function can draw multiple ellipsoids with a single patch command.
% Significant speedup can be achieved by a single patch compared to ellipse
% per patch.
%
% handle = SWPLOT.ELLIPSOID(hAxis,...
%
% plots to the selected axis
%
% Input:
%
% hAxis     Axis handle.
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
    % first input figure handle
    hAxis = varargin{1};
    R0  = varargin{2};
    T    = varargin{3};
    if nargin > 3
        mesh = varargin{4};
    else
        mesh = [];
    end
else
    hAxis = gca;
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

% plot multiple ellipse
nEllipse = size(R0,2);

if size(R0,1)~=3 || any(size(T)~=[3 3 nEllipse])
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

% create patch
handle = patch(hAxis,'Vertices',V,'Faces',F,'FaceLighting','flat',...
    'EdgeColor','none','FaceColor','r','Tag','ellipsoid');

end