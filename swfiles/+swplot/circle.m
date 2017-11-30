function hPatch = circle(varargin)
% creates a 3D circle surface patch
% 
% ### Syntax
% 
% `hPatch = swplot.circle(r0, n, R)`
% 
% `hPatch = swplot.circle(r0, n, R, nPatch)`
%
% `hPatch = swplot.circle(handle, ...)`
%
% ### Description
% 
% `hPatch = swplot.circle(r0, n, R)` creates a triangulated patch of a
% surface of a circle in 3D, defined by the center position, normal vector
% and radius.
%  
% `hPatch = swplot.circle(handle, ...)` adds the patch object to a given axis
% if `handle` is an axis handle or adds the arrow to an existing
% [matlab.patch] object, if the given `handle` points to a patch object.
%  
% 
% ### Examples
%
% Draw 100 random unit circle surfaces with center at $(0,0,0)$ and random
% normal vector.
%
% ```
% >>swplot.figure
% >>N = 100
% >>swplot.circle(zeros(3,N),2*rand(3,N)-1,1)
% >>swplot.zoom(30)% 
% >>snapnow
% ```
%
% ### Input Arguments
% 
% `handle`
% : Handle of an axis or triangulated patch object. In case of patch
%   object, the constructed faces will be added to the existing object.
% 
% `r0`
% : Center position of the circle in a column vector. Multiple circles can
%   be defined using a matrix with dimensions of $[3\times n_{obj}]$ where
%   each column defines a circle center.
% 
% `n`
% : Column vector with 3 elements, normal to the circle surface. Multiple
%   circles can be defined using a matrix with the same dimensions as `r0`
%   parameter.
% 
% `R`
% : Radius of the circle, scalar or row vector with $n_{obj}$ number of
%   elements.
% 
% `nPatch`
% : Number of points on the circle circumference, default value is stored in
%   `swpref.getpref('npatch')`. The generated patch will contain
%   $n_{patch}$ number of faces and vertices.
% 
% ### See Also
% 
% [swplot.cylinder]
%

if nargin == 0
    help swplot.circle
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
    r0      = varargin{2};
    n       = varargin{3};
    R       = varargin{4};
    nArgExt = nargin-4;
    argExt  = {varargin{5:end}};
    
else
    hAxis   = gca;
    hPatch  = [];
    r0      = varargin{1};
    n       = varargin{2};
    R       = varargin{3};
    nArgExt = nargin-3;
    argExt  = {varargin{4:end}};
end

pref = swpref;

if nArgExt > 0
    N = argExt{1};
else
    N = pref.npatch;
end

if numel(r0) == 3
    r0 = r0(:);
    n  = n(:);
end

nObject = size(r0,2);

% normal vectors to the cylinder axis
a = cross(n,repmat([0;0;1],[1 nObject]));

% index of zero normal vectors
zIdx = find(sum(abs(a),1)==0);
% try another normal vector for these
if ~isempty(zIdx)
    a(:,zIdx) = cross(n(:,zIdx),repmat([0;1;0],[1 numel(zIdx)]));
end

b = cross(n,a);
a = bsxfun(@rdivide,a,sqrt(sum(a.^2,1)));
b = bsxfun(@rdivide,b,sqrt(sum(b.^2,1)));

phi    = permute(linspace(0,2*pi,N+1),[1 3 2]);
phi    = phi(1,1,1:N);
cPoint = R*(bsxfun(@times,a,cos(phi))+bsxfun(@times,b,sin(phi)));

% vertices
V = reshape(permute(bsxfun(@plus,cPoint,r0),[1 3 2]),3,[])';

% faces
L = (2:(N-1))';
F = [ones(N-2,1) L mod(L,N)+1];

F = reshape(permute(bsxfun(@plus,F,permute((0:(nObject-1))*N,[1 3 2])),[1 3 2]),[],3);

% color data
C = repmat([1 0 0],[size(F,1) 1]);
% default transparency
A = ones(size(F,1),1);

if isempty(hPatch)
    % create patch
    hPatch = patch('Parent',hAxis,'Vertices',V,'Faces',F,'FaceLighting','flat',...
        'EdgeColor','none','FaceColor','flat','Tag','circle','FaceVertexCData',C,...
        'AlphaDataMapping','none','FaceAlpha','flat','FaceVertexAlphaData',A);
else
    % add to existing patch
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
    hPatch = repmat(hPatch,[1 nObject]);
end

end