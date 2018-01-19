function hPatch = ellipsoid(varargin)
% creates a 3D ellipsoid patch
% 
% ### Syntax
% 
% `hPatch = swplot.ellipsoid(R0, T)`
%
% `hPatch = swplot.ellipsoid(R0, T, nMesh)`
% 
% `hPatch = swplot.ellipsoid(handle, ...)`
%
% ### Description
% 
% `hPatch = swplot.ellipsoid(R0,T)` creates multiple ellipsoids with a
% single [matlab.patch] command. The ellipsoids are defined by the position
% of the center and a $[3\times 3]$ matrix, a qudratic form.
%
% Significant speedup can be achieved using a single patch command to
% generate many ellipsoids compared to drawing single ellipse per patch.
%  
% `hPatch = swplot.ellipsoid(R0, T, nMesh)` defines the size of the mesh
% that defines the surface.
%
% `hPatch = swplot.ellipsoid(handle, ...)` adds the generated patch object
% to a given axis if `handle` is an axis handle or adds the ellipsoids to
% an existing [matlab.patch] object, if the given `handle` points to a
% patch object.
% 
% ### Input Arguments
% 
% `handle`
% : Handle of an axis or triangulated patch object. In case of patch
%   object, the constructed faces will be added to the existing object.
% 
% `R0`
% : Center of the ellipsoids stored in a column vector with 3 elements or a
%   matrix with dimensions of $[3\times n_{obj}]$ when multiple ellipsoids
%   are defined at once.
% 
% `T`
% : Transformation matrix that transforms a unit sphere to the desired
%   ellipsoid by applying: `R' = T(:,:,i)*R`. In case of multiple
%   ellipsoids the parameter is stored in a matrix with dimensions of
%   $[3\times 3\times n_{obj}]$.
% 
% `nMesh`
% : Mesh of the ellipse surface, a triangulation class object or an
%   integer that used to generate an icosahedron mesh with $n_{mesh}$
%   number of additional subdivision into triangles. Default value is stored in
%   `swpref.getpref('nmesh')`, see also [swplot.icomesh].
% 
% ### See Also
% 
% [matlab.triangulation] \| [swplot.icomesh]
%
pref = swpref;

if nargin == 0
    swhelp swplot.ellipsoid
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
    mesh = min(mesh,pref.maxmesh);
    % generate mesh
    mesh = swplot.icomesh(mesh);
end

if numel(R0) == 3
    R0 = R0(:);
end

% plot multiple ellipse
nObject = size(R0,2);

sT = [size(T) 1];
if size(R0,1)~=3 || any(sT(1:3)~=[3 3 nObject])
    error('ellipsoid:WrongInput','Matrices have incompatible dimensions!')
end

if isa(mesh,'TriRep')
    V0 = mesh.X;
    F0 = mesh.Triangulation;
elseif isa(mesh,'triangulation')
    V0 = mesh.Points;
    F0 = mesh.ConnectivityList;
else
    error('ellipsoid:WrongInput','The given data is not a Matlab triangulation/TriRep class!')
end

NV = size(V0,1);

% fast vertices
%V = reshape(permute(bsxfun(@plus,sum(bsxfun(@times,T,permute(V0,[3 2 4 1])),1),permute(R0,[3 1 2])),[4 3 2 1]),[],3);
V = reshape(permute(bsxfun(@plus,sum(bsxfun(@times,T,permute(V0,[3 2 4 1])),2),permute(R0,[1 3 2])),[4 3 1 2]),[],3);
% fast faces
F = reshape(bsxfun(@plus,permute(F0,[1 3 2]),((1:nObject)-1)*NV),[],3);

% default red color
C = repmat([1 0 0],[size(F,1) 1]);
% default transparency
A = ones(size(F,1),1);

if isempty(hPatch)
    % create patch
    hPatch = patch('Parent',hAxis,'Vertices',V,'Faces',F,'FaceLighting','flat',...
        'EdgeColor','none','FaceColor','flat','Tag','ellipsoid','AlphaDataMapping','none',...
        'FaceAlpha','flat','FaceVertexAlphaData',A,'FaceVertexCData',C);
else
    V0 = get(hPatch,'Vertices');
    F0 = get(hPatch,'Faces');
    C0 = get(hPatch,'FaceVertexCData');
    A0 = get(hPatch,'FaceVertexAlphaData');
    % number of existing faces
    nV0 = size(V0,1);
    set(hPatch,'Vertices',[V0;V],'Faces',[F0;F+nV0],'FaceVertexCData',[C0;C],...
        'FaceVertexAlphaData',[A0;A]);
end

if strcmp(get(hAxis,'Tag'),'swaxis')
    % replicate the arrow handle to give the right number of added objects
    % for swplot figure
    hPatch = repmat(hPatch,[1 nObject]);
end

end