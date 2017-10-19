function [tuv, fIdx] = raytriangle(V,F,ray)
% finds if a ray crosses a triangle
% 
% ### Syntax
% 
% `swplot.raytriangle(V,F,ray)`
% 
% ### Description
% 
% `swplot.raytriangle(V,F,ray)` determines if a given ray crosses within a
% triangle defined by its vertex coordinates. The code is optimised for
% a single ray and multiple triangles.
% 
% ### Input Arguments
% 
% `V`
% : Vertex positions in a matrix with dimensions $[n_{vertex}\times 3]$.
% 
% `F`
% : Faces in a matrix with dimensions $[n_{face}\times 3]$, where 
%       $max(F) = n_{vertex}$.
% 
% `ray`
% : Definition of the ray via 2 points in space, stored in a matrix with
%   dimensions of $[2\times 3]$. The two points of the ray, $P_1$ and
%   $P_2$, are stored in the first and second rows respectively. The ray
%   points from $P_1$ to $P_2$.
%

% number of faces
nF = size(F,1);

% subtract the origin of the triangles
E1 = V(F(:,2),:)-V(F(:,1),:);
E2 = V(F(:,3),:)-V(F(:,1),:);

% define line with origin and unit vector
O = ray(1,:);
D = ray(2,:)-O;
D = D/norm(D);

P = cross(repmat(D,[nF 1]),E2);
% det defines whether the ray lies in the plane of triangle
det = sum(E1.*P,2);

T = bsxfun(@minus,O,V(F(:,1),:));

% u coordinate
u = sum(T.*P,2);

% test u-coordinate
inIdxU = u>=0 & u<=det;

% continue calculation only on faces that are inside
Q = cross(T(inIdxU,:),E1(inIdxU,:));
v = sum(repmat(D,[size(Q,1),1]).*Q,2);

% test v-coordinate
detv = det(inIdxU);
inIdxV = v>=0 & v<=detv;

E2u = E2(inIdxU,:);
t   = sum(E2u(inIdxV,:).*Q(inIdxV,:),2);

invDet = 1./detv(inIdxV);

uu = u(inIdxU);

t = t.*invDet;
u = uu(inIdxV).*invDet;
v = v(inIdxV).*invDet;

% find the intersecting faces
fIdx = find(inIdxU);
fIdx = fIdx(inIdxV);

% give tuv for intersecting rays
tuv = [t u v];

end