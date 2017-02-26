function hPatch = orbital(varargin)
% draw electron orbitals
%
% hPatch = SWPLOT.ORBITAL(qLabel,R0,T,scale,nPatch)
%
% The function can draw electron orbitals with a single patch command.
%
% hPatch = SWPLOT.ORBITAL(handle,...
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
% qLabel    Label of the orbital in a string, corresponding quantum
%           numbers:
%                                   n l m pm
%               's'         qNum = [1 0 0  0]
%               'p_x'       qNum = [2 1 1  1]
%               'p_y'       qNum = [2 1 0 -1]
%               'p_z'       qNum = [2 1 0  0]
%               'd_xy'      qNum = [3 2 2 -1]
%               'd_xz'      qNum = [3 2 1  1]
%               'd_yz'      qNum = [3 2 1 -1]
%               'd_z2'      qNum = [3 2 0  0]
%               'd_x2-y2'   qNum = [3 2 2  1]
%           where the quantum numbers are in a vector (n, l, m, {pm}),
%           where pm defines the optional linear combination of the +m and
%           -m orbitals (PSI is the wave function):
%               PSI = PSI(n,l,m) + pm*PSI(n,l,-m)
%           If pm is +1,-1 or m=0 the wave fuction is real, otherwise
%           complex.
% R0        Center of the orbital stored in a matrix with dimensions of
%           [3 nOrbital].
% T         Rotation matrix that rotates the orbital via:
%           O' = T(:,:,i)*O*T(:,:,i)'
%           Dimensions are [3 3 nOrbital].
% scale     Scale of the orbital. Default scale is the true Hydgroden
%           orbital size in Angstrom. Default value is 1.
% mesh      Mesh of the ellipse surface, a triangulation class object or an
%           integer that used to generate an icosahedron mesh with #mesh
%           number of additional triangulation. Default value is stored in
%           swpref.getpref('nmesh')
%
% See also TRIANGULATION, SWPLOT.ICOMESH.
%

if nargin == 0
    help swplot.orbital
    return
end

if numel(varargin{1}) == 1 && ~ischar(varargin{1})
    % first input figure/patch handle
    if strcmp(get(varargin{1},'Type'),'axes')
        hAxis  = varargin{1};
        hPatch = [];
    else
        hAxis  = gca;
        hPatch = varargin{1};
    end
    qLabel = varargin{2};
    R0     = varargin{3};
    T      = varargin{4};
    scale  = varargin{5};
    if nargin > 5
        nPatch = varargin{6};
    else
        nPatch = [];
    end
    
else
    hAxis  = gca;
    hPatch = [];
    qLabel = varargin{1};
    R0     = varargin{2};
    T      = varargin{3};
    scale  = varargin{4};
    if nargin > 4
        nPatch = varargin{5};
    else
        nPatch = [];
    end
end

if isempty(nPatch)
    nPatch = 20;
end

% generate mesh of the orbital
[mesh, S] = swplot.orbmesh(qLabel,'nPatch',nPatch);

if numel(R0) == 3
    R0 = R0(:);
end

% plot multiple orbitals
nObject = size(R0,2);

sT = [size(T) 1];
if size(R0,1)~=3 || any(sT(1:3)~=[3 3 nObject])
    error('orbital:WrongInput','Matrices have incompatible dimensions!')
end

if isa(mesh,'TriRep')
    V0 = mesh.X;
    F0 = mesh.Triangulation;
elseif isa(mesh,'triangulation')
    V0 = mesh.Points;
    F0 = mesh.ConnectivityList;
else
    error('orbital:WrongInput','The given data is not a Matlab triangulation/TriRep class!')
end

% scale the orbital
V0 = V0*scale;

NV = size(V0,1);

% fast vertices
V = reshape(permute(bsxfun(@plus,sum(bsxfun(@times,T,permute(V0,[3 2 4 1])),2)...
    ,permute(R0,[1 3 2])),[4 3 1 2]),[],3);
% fast faces
F = reshape(bsxfun(@plus,permute(F0,[1 3 2]),((1:nObject)-1)*NV),[],3);

% default red-blue color
C = [1 0 0;0 0 1];
C = repmat(C(S,:),nObject,1);
% default transparency
A = ones(size(F,1),1);

if isempty(hPatch)
    % create patch
    hPatch = patch('Parent',hAxis,'Vertices',V,'Faces',F,'FaceLighting','flat',...
        'EdgeColor','none','FaceColor','flat','Tag','orbital','AlphaDataMapping','none',...
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