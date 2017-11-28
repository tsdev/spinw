function hPatch = arrow(varargin)
% creates a 3D arrow patch
% 
% ### Syntax
% 
% `hPatch = swplot.arrow(rStart, rEnd, R, alpha, lHead)`
%
% `hPatch = swplot.arrow(rStart, rEnd, R, alpha, lHead, nPatch)`
% 
% `hPatch = swplot.arrow(handle, ...)`
%
% ### Description
% 
% `hPatch = swplot.arrow(rStart, rEnd, R, alpha, lHead)` draws 3D arrows
% between a given start and end position. The arrows will be a triangulated
% [matlab.patch] object.
%  
% `hPatch = swplot.arrow(rStart, rEnd, R, alpha, lHead, nPatch)` creates 
% arrows with $5 n_{patch}$ number of patch faces per arrow.
%
% `hPatch = swplot.arrow(handle, ...)` adds the generated patch object to a
% given axis if `handle` is an axis handle or adds the arrows to an
% existing [matlab.patch] object, if the given `handle` points to a patch
% object.
%  
% ### Examples
%
% Draw a 100 random arrows in the $(\pm 1,\pm 1,\pm 1)$ cube:
%
% ```
% >>swplot.figure
% >>N = 100
% >>swplot.arrow(2*rand(3,N)-1,2*rand(3,N)-1,0.01,30,0.05)
% >>swplot.zoom(40)
% >>snapnow
% ```
%
% ### Input Arguments
% 
% `handle`
% : Handle of an axis or triangulated patch object. In case of patch
%   object, the constructed faces will be added to the existing object.
% 
% `rStart`
% : Coordinates of the arrow starting point, one vector per arrow in a
%   matrix with dimensions of $[3\times n_{obj}]$.
% 
% `rEnd`
% : Coordinates of the arrow end point, one vector per arrow in a
%   matrix with dimensions of $[3\times n_{obj}]$.
% 
% `R`
% : Radius of the arrow body, scalar.
% 
% `alpha`
% : Angle of the head in degree.
% 
% `lHead`
% : Length of the head.
% 
% `nPatch`
% : Number of points on the circle of the body, default value is stored in
%   `swpref.getpref('npatch')`. The final patch object will have
%   $5n_{patch}$ number of faces and $3n_{patch}$ number of vertices.
% 
% ### See Also
% 
% [swplot.cylinder]
%

if nargin == 0
    help swplot.arrow
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
    rStart  = varargin{2};
    rEnd    = varargin{3};
    R       = varargin{4};
    alpha   = varargin{5};
    lHead   = varargin{6};
    if nargin > 6
        nPatch = varargin{7};
    else
        nPatch = [];
    end
else
    hAxis   = gca;
    hPatch  = [];
    rStart  = varargin{1};
    rEnd    = varargin{2};
    R       = varargin{3};
    alpha   = varargin{4};
    lHead   = varargin{5};
    if nargin > 5
        nPatch = varargin{6};
    else
        nPatch = [];
    end
    
end

if isempty(nPatch)
    nPatch = swpref.getpref('npatch',[]);
end

if numel(rStart)==3
    rStart = rStart(:);
    rEnd   = rEnd(:);
end

rArrow  = rEnd - rStart;
endHead = rStart+bsxfun(@times,1-lHead./sqrt(sum(rArrow.^2,1)),rArrow);
rHead   = lHead * tand(alpha);

% number of cylinder segments
nObject = size(rArrow,2);

% normal vectors to the cylinder axis
nArrow1 = cross(rArrow,repmat([0;0;1],[1 nObject]));

% index of zero normal vectors
zIdx = find(sum(abs(nArrow1),1)==0);
% try another normal vector for these
if ~isempty(zIdx)
    nArrow1(:,zIdx) = cross(rArrow(:,zIdx),repmat([0;1;0],[1 numel(zIdx)]));
end

nArrow2 = cross(rArrow,nArrow1);
nArrow1 = bsxfun(@rdivide,nArrow1,sqrt(sum(nArrow1.^2,1)));
nArrow2 = bsxfun(@rdivide,nArrow2,sqrt(sum(nArrow2.^2,1)));

phi     = permute(linspace(0,2*pi,nPatch+1),[1 3 2]);
phi     = phi(1,1,1:nPatch);
cPoint  = R*(bsxfun(@times,nArrow1,cos(phi))+bsxfun(@times,nArrow2,sin(phi)));
cPoint1 = bsxfun(@plus,cPoint,rStart);
cPoint2 = bsxfun(@plus,cPoint,endHead);
cPoint3 = bsxfun(@plus,cPoint*rHead/R,endHead);

% vertices
V = reshape(permute(cat(3,cPoint1,cPoint2,cPoint3,rEnd),[1 3 2]),3,[])';

% faces
% caps
L  = (2:(nPatch-1))';
F1 = [ones(nPatch-2,1) L mod(L,nPatch)+1];
% body
L  = (1:nPatch)';
F2 = [L L+nPatch mod(L,nPatch)+nPatch+1 L mod(L,nPatch)+1 mod(L,nPatch)+nPatch+1];
F2 = reshape(F2',3,[])';
% head back
F3 = F1+2*nPatch;
% head
F4 = [L mod(L,nPatch)+1 repmat(nPatch+1,[nPatch 1])] + 2*nPatch;

F = [F1;F2;F3;F4];

F = reshape(permute(bsxfun(@plus,F,permute((0:(nObject-1))*(3*nPatch+1),[1 3 2])),[1 3 2]),[],3);

% color data
C = repmat([1 0 0],[size(F,1) 1]);
% default transparency
A = ones(size(F,1),1);

if isempty(hPatch)
    % create patch
    hPatch = patch('Parent',hAxis,'Vertices',V,'Faces',F,'FaceLighting','flat',...
        'EdgeColor','none','FaceColor','flat','Tag','arrow','AlphaDataMapping','none',...
        'FaceAlpha','flat','FaceVertexAlphaData',A,'FaceVertexCData',C);
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