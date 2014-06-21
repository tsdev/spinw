function [handle] = sw_cone(rBase, rTop, R, N)
% creates the surface of a cone in 3 dimension
%
% [handle] = SW_CONE(rBase, rTop, R, N)
%
% Input:
%
% rBase    Coordinates of the middle of the base, dimensions are [3 1].
% rTop     Coordinates of the top, dimensions are [3 1].
% R        Radius of the base of the cone.
% N        number of points of the surface mesh
%
% See also SW_CIRCLE, SW_CIRCLESURF, SW_CYLINDER, SW_ARROW.
%

if nargin == 0
    help sw_cone;
    return
end

n = rTop - rBase;

if any(cross(n,[0; 0; 1]))
    a = cross(n,[0; 0; 1]);
else
    a = cross(n,[0; 1; 0]);
end

b = cross(n,a);
a = a/norm(a);
b = b/norm(b);

phi = linspace(0,2*pi,N);

points = R*(a*cos(phi)+b*sin(phi))+repmat(rBase,1,N);

X = [points(1,:); rTop(1)*ones(1,N)];
Y = [points(2,:); rTop(2)*ones(1,N)];
Z = [points(3,:); rTop(3)*ones(1,N)];

handle = surface(X,Y,Z,'LineStyle','none');

end