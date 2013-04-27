function [handle] = sw_circlesurf(r0, n, R, N)
% [handle] = SW_CIRCLESURF(r0, n, R, N) creates the filled of a circle in 
% 3 dimensions.
%
% Input:
%
% r0    Center of the circle, dimensions are [3 1].
% n     Vector normal to the circle surface, dimensions are [3 1].
% R     Radius of the circle.
% N     Number of points on the curve.
%
% See also SW_CIRCLE, SW_CONE, SW_CYLINDER.
%

if any(cross(n,[0; 0; 1]))
    a = cross(n,[0; 0; 1]);
else
    a = cross(n,[0; 1; 0]);
end

b = cross(n,a);
a = a/norm(a);
b = b/norm(b);

phi = linspace(0,2*pi,N);

edge = R*(a*cos(phi)+b*sin(phi))+repmat(r0,1,N);


X = [edge(1,:); r0(1)*ones(1,N)];
Y = [edge(2,:); r0(2)*ones(1,N)];
Z = [edge(3,:); r0(3)*ones(1,N)];

handle = surface(X,Y,Z,'LineStyle','none');

end