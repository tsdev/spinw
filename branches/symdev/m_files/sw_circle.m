function points = sw_circle(r0, n, R, N)
% points = SW_CIRCLE(r0, n, R, N) creates array containing the 3D points of
% the circle line.
% r0    Center of circle, dimensions are [3 1].
% n     Normal to the circle surface, dimensions are [3 1].
% R     Radius of the circle.
% N     Number of points on the curve.
%
% See also SW_CONE, SW_CIRCLESURF, SW_CYLINDER.
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

points = R*(a*cos(phi)+b*sin(phi))+repmat(r0,1,N);

end