function [handle] = sw_arrow(rStart, rEnd, R, alpha, lHead, N)
% draws a 3D arrow
%
% [handle] = SW_ARROW(rStart, rEnd, R, alpha, lHead, N)
%
% Input:
%
% rStart  Coordinate of the starting point.
% rEnd    Coordinate of the end point.
% R       Radius of the arrow body.
% alpha   Angle of the head in degrees.
% lHead   Length of the head.
% N       Number of the points of the surface mesh.
%

if nargin == 0
    help sw_arrow
    return
end

rStart = rStart(:);
rEnd   = rEnd(:);

alpha    = alpha *pi/180;
end_head = rStart+(1-lHead/norm(rEnd-rStart))*(rEnd-rStart);
R_head   = lHead * tan(alpha);

%body
handle = sw_cylinder(rStart,end_head+1e-2*(rEnd-rStart),R,N,0);
%annulus for head
handle(end+1) = sw_circlesurf(end_head, rEnd-rStart, R_head, N);
%cone for head
handle(end+1) = sw_cone(end_head, rEnd, R_head, N);
%closing the begin of body
handle(end+1) = sw_circlesurf(rStart, rEnd-rStart, R, N);

end