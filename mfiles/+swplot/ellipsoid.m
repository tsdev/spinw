function hEllipse = ellipsoid(R0,T,mesh, varargin)
% draw ellipsoid
%
% hEllipse = SWPLOT.ELLIPSOID(R0,T,mesh, 'Option1', value1, ...)
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
% Options:
%
% Any property of the patch class can be set.
%
% See also TRIANGULATION, SWPLOT.ICOMESH.
%

if isempty(mesh)
    mesh = 1;
end

if isnumeric(mesh)
    % generate mesh
    mesh = swplot.icomesh(mesh);
end

% ellipse points
X = bsxfun(@plus,T*mesh.Points',R0(:))';

% create patch
hEllipse = trimesh(mesh.ConnectivityList,X(:,1),X(:,2),X(:,3),'FaceLighting','flat','EdgeAlpha',0,varargin{:});

end