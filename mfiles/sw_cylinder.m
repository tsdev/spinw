function [handle] = sw_cylinder(r1,r2,R,N, dash)
% creates the surface of a cylinder in 3D
%
% [handle] = SW_CYLINDER(r1,r2,R,N, dash) 
%
% Input:
%
% r1    Coordinates of the starting point, dimensions are [3 1].
% r2    Coordinates of the end point, dimensions are [3 1].
% R     Radius of the cylinder.
% N     Number of points of the surface mesh.
% dash  If defined either as a single number or a vector [dash N1 N2], a
%       dashed cylinder surface will be drawn, with gaps as follows:
%           [N1*R(gap)*dash N2*R(dash)*dash].
%       Default for dash, N1 and N2 are 0, 4 and 1 respectively.
%
% See also SW_CIRCLE, SW_CONE, SW_CIRCLESURF.
%

if nargin == 0
    help sw_cylinder;
    return
end

if nargin < 5
    dash = 0;
    N1 = 4;
    N2 = 1;
elseif numel(dash) == 1
    N1 = 4;
    N2 = 1;
else
    N1 = dash(2);
    N2 = dash(3);
    dash = dash(1);
end
    
if dash == 0
    nDash = 0;
else
    nDash = floor(norm(r2-r1)/((N1+N2)*R*dash));
end

% axis of the cylinder
ax = (r2-r1)/norm(r2-r1);

handle = zeros(1,nDash);

rf = r1;

for ii = 1:nDash
    ri = ax*((ii-1)*(N1+N2)*R*dash+N1*R*dash)+r1;
    rf = ri + N2*ax*R*dash;
    
    c1 = sw_circle(ri, ax, R, N);
    c2 = sw_circle(rf, ax, R, N);
    
    X = [c1(1,:); c2(1,:)];
    Y = [c1(2,:); c2(2,:)];
    Z = [c1(3,:); c2(3,:)];
    
    handle(ii) = surface(X,Y,Z,'LineStyle','none');
    
end

if norm(r2-rf)>(2*R*dash)
    if nDash > 0
        ri = ax*(nDash*(N1+N2)*R*dash+N1*R*dash)+r1;
    else
        ri = r1;
    end
    
    c1 = sw_circle(ri, ax, R, N);
    c2 = sw_circle(r2, ax, R, N);
    
    X = [c1(1,:); c2(1,:)];
    Y = [c1(2,:); c2(2,:)];
    Z = [c1(3,:); c2(3,:)];
    
    handle(nDash+1) = surface(X,Y,Z,'LineStyle','none');
    
end


end