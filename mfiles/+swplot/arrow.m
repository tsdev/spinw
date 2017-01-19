function hPatch = arrow(varargin)
% draws a 3D arrow using patch
%
% hPatch = SWPLOT.ARROW(rStart, rEnd, R, alpha, lHead, {nPatch})
%
% hPatch = SWPLOT.ARROW(handle,...)
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
% rStart    Coordinate of the starting point.
% rEnd      Coordinate of the end point.
% R         Radius of the arrow body.
% alpha     Angle of the head in degrees.
% lHead     Length of the head.
% nPatch    Number of points on the curve, default value is stored in
%           swpref.getpref('npatch').
%
% See also SWPLOT.CYLINDER.
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