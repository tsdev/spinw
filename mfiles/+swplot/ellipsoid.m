function hEllipse = ellipsoid(R0,T,mesh)
% draw ellipsoid
%
% hEllipse = SWPLOT.ELLIPSOID(R0,T,mesh)
%
% Input:
%
% R0        Center of the ellipsoid stored in a 3 element vector.
% T         Transformation matrix that transforms a unit sphere to the
%           ellipse via: R' = T*R
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

if isempty(mesh)
    mesh = 1;
end

if isnumeric(mesh)
    % limit the largest mesh to plot to avoid slowing down Matlab too much
    mesh = min(mesh,swpref.getpref('maxmesh',[]));
    % generate mesh
    mesh = swplot.icomesh(mesh);
end

% ellipse points
X = bsxfun(@plus,T*mesh.Points',R0(:))';

% create patch
hEllipse = trimesh(mesh.ConnectivityList,X(:,1),X(:,2),X(:,3),'FaceLighting','flat','EdgeColor','none','FaceColor','r');

end