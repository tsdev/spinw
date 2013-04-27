function [handle] = sw_cylinder(r1,r2,R,N)
% [handle] = SW_CYLINDER(r1,r2,R,N) creates the surface of the cylinder.
%
% Input:
%
% r1    Coordinates of the starting point, dimensions are [3 1].
% r2    Coordinates of the end point, dimensions are [3 1].
% R     Radius of the cylinder.
% N     Number of points of the surface mesh.
%
% See also SW_CIRCLE, SW_CONE, SW_CIRCLESURF.
%

c1 = sw_circle(r1, r2-r1, R, N);
c2 = sw_circle(r2, r2-r1, R, N);

X = [c1(1,:); c2(1,:)];
Y = [c1(2,:); c2(2,:)];
Z = [c1(3,:); c2(3,:)];

handle = surface(X,Y,Z,'LineStyle','none');

end