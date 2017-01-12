function hPatch = circle(varargin)
% creates a circle surface in 3 dimensions
%
% hPatch = SWPLOIT.CIRCLE(r0, n, R, {N})
%
% hPatch = SWPLOIT.CIRCLE(handle,...)
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
% r0        Center of the circle, vector with three elements.
% n         Vector normal to the circle surface, vector with three elements.
% R         Radius of the circle.
% N         Number of points on the curve, default value is stored in
%           swpref.getpref('npatch').
%
% Example:
%
% swplot.circle(zeros(3),eye(3),1,100)
%
% See also SWPLOT.CYLINDER.
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

if nArgExt > 0
    N = argExt{1};
else
    N = swpref.getpref('npatch',[]);
end

if numel(r0) == 3
    r0 = r0(:);
    n  = n(:);
end

nCircle = size(r0,2);

% normal vectors to the cylinder axis
a = cross(n,repmat([0;0;1],[1 nCircle]));

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

F = reshape(permute(bsxfun(@plus,F,permute((0:(nCircle-1))*N,[1 3 2])),[1 3 2]),[],3);

% color data
C = repmat([1 0 0],[size(F,1) 1]);

if strcmp(get(hAxis,'Tag'),'swaxis')
    % make sure we are on the plot axis of an swobject
    % add object to the existing triangular patch
    hFigure = get(hAxis,'Parent');
    hPatch = getappdata(hFigure,'facepatch');
end

if isempty(hPatch)
    % create patch
    hPatch = patch(hAxis,'Vertices',V,'Faces',F,'FaceLighting','flat',...
        'EdgeColor','none','FaceColor','flat','Tag','circle','FaceVertexCData',C);
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