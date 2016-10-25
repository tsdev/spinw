function points = sw_ellipse(r0, n, R, N, v1, v2)
% creates an array of the 3D coordinates of the ellipse circumference
%
% points = SW_ELLIPSE(r0, n, R, N, {v1, v2})
%
% Input:
%
% r0    Center of ellipse, dimensions are [3 1].
% n     Normal vector to the ellipse surface, dimensions are [3 1].
% R     Minor and major radius of the ellipse, dimensions are [2 1].
% N     Number of points on the curve.
% v1,v2 Optional parameter, two vectors pointing parallel to the minor and
%       major radius of the ellipse, dimensions are [3 1].
%
% See also SW_CONE, SW_CIRCLE, SW_CIRCLESURF, SW_CYLINDER.
%

if nargin == 0
    help sw_ellipse;
    return
end

if nargin < 6
    if any(cross(n,[0; 0; 1]))
        v1 = cross(n,[0; 0; 1]);
    else
        v1 = cross(n,[0; 1; 0]);
    end
    
    v2 = cross(n,v1);
    v1 = v1/norm(v1);
    v2 = v2/norm(v2);
end

phi = linspace(0,2*pi,N);

points = R(1)*v1*cos(phi) + R(2)*v2*sin(phi) + repmat(r0,1,N);

end