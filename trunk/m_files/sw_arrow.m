function [handle] = sw_arrow(rStart, rEnd, R, alpha, lHead, N)
% [handle] = SW_ARROW(rStart, rEnd, R, alpha, lHead, N) plots a 3D arrow.
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

rStart = rStart(:);
rEnd   = rEnd(:);

alpha    = alpha *pi/180;
end_head = rStart+(1-lHead/norm(rEnd-rStart))*(rEnd-rStart);
R_head   = lHead * tan(alpha);
handle   = zeros(4,1);

%body
handle(1) = sw_cylinder(rStart,end_head+1e-2*(rEnd-rStart),R,N);
%annulus for head
handle(2) = sw_circlesurf(end_head, rEnd-rStart, R_head, N);
%cone for head
handle(3) = sw_cone(end_head, rEnd, R_head, N);
%closing the begin of body
handle(4) = sw_circlesurf(rStart, rEnd-rStart, R, N);

end