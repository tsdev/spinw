function points = sw_circle(r0, n, R, N)
% creates an array of the 3D coordinates of the circle circumference
%
% points = SW_CIRCLE(r0, n, R, N) 
%
% r0    Center of circle, dimensions are [3 1].
% n     Normal to the circle surface, dimensions are [3 1].
% R     Radius of the circle.
% N     Number of points on the curve.
%
% See also SW_CONE, SW_CIRCLESURF, SW_CYLINDER, SW_ARROW.
%

% $Name: SpinW$ ($Version: 2.1$)
% $Author: S. Toth$ ($Contact: sandor.toth@psi.ch$)
% $Revision: 238 $ ($Date: 07-Feb-2015 $)
% $License: GNU GENERAL PUBLIC LICENSE$

if nargin == 0
    help sw_circle;
    return
end

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