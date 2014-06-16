function [handle] = sw_ellcylinder(r1,r2,R,N,v1,v2, cap)
% [handle] = SW_ELLCYLINDER(r1,r2,R,N,{v1, v2, cap}) creates the surface of
% an elliptic cylinder.
%
% Input:
%
% r1    Coordinates of the center of the base ellipse, dimensions are [3 1].
% r2    Coordinates of the center of the top cylinder, dimensions are [3 1].
% R     Minor and major radius of the ellipse, dimensions are [2 1].
% N     Number of points of the surface mesh.
% v1,v2 Optional nomal vectors pointing along the minor and major radius
%       vector, dimensions are [3 1].
% cap   Optional input, if true the elliptical cylinder is closed with
%       caps. Default is false.
%
% See also SW_CIRCLE, SW_ELLIPSE, SW_CONE, SW_CIRCLESURF, SW_CYLINDER.
%

if nargin == 0
    help sw_ellcylinder;
    return
end

% By default no cap is drawn
if nargin < 7
    cap = false;
end

% axis of the elliptical cylinder
ax = (r2-r1)/norm(r2-r1);

if nargin < 6
    c1 = sw_ellipse(r1, ax, R, N);
    c2 = sw_ellipse(r2, ax, R, N);
else
    c1 = sw_ellipse(r1, ax, R, N, v1, v2);
    c2 = sw_ellipse(r2, ax, R, N, v1, v2);
end

X = [c1(1,:); c2(1,:)];
Y = [c1(2,:); c2(2,:)];
Z = [c1(3,:); c2(3,:)];

handle(1) = surface(X,Y,Z);

if cap
    X = [c1(1,:); r1(1)*ones(1,N)];
    Y = [c1(2,:); r1(2)*ones(1,N)];
    Z = [c1(3,:); r1(3)*ones(1,N)];
    handle(2) = surface(X,Y,Z);
    
    X = [c2(1,:); r2(1)*ones(1,N)];
    Y = [c2(2,:); r2(2)*ones(1,N)];
    Z = [c2(3,:); r2(3)*ones(1,N)];
    handle(3) = surface(X,Y,Z);
    
end

set(handle,'LineStyle','none');

end