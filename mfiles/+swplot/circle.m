function [handle] = circle(r0, n, R, N)
% creates a circle surface in 3 dimensions
%
% handle = SWPLOIT.CIRCLE(r0, n, R, {N})
%
% Input:
%
% r0        Center of the circle, vector with three elements.
% n         Vector normal to the circle surface, vector with three elements.
% R         Radius of the circle.
% N         Number of points on the curve, default value is stored in 
%           swpref.getpref('npatch').
%
% See also SWPLOT.CYLINDER.
%

if nargin == 0
    help swplot.circle
    return
end

if nargin < 4
    N = swpref.getpref('npatch',[]);
end

r0 = repmat(r0(:),1,N);
n  = n(:);

if any(cross(n,[0; 0; 1]))
    a = cross(n,[0; 0; 1]);
else
    a = cross(n,[0; 1; 0]);
end

b = cross(n,a);
a = a/norm(a);
b = b/norm(b);

phi = linspace(0,2*pi,N);

edge = mat2cell(R*(a*cos(phi)+b*sin(phi))+r0,ones(1,3),N);

handle = patch(edge{:},'FaceLighting','flat','EdgeColor','none','FaceColor','r');

end