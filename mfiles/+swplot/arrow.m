function handle = arrow(rStart, rEnd, R, alpha, lHead, N)
% draws a 3D arrow
%
% [handle] = SWPLOT.ARROW(rStart, rEnd, R, alpha, lHead, {N})
%
% Input:
%
% rStart    Coordinate of the starting point.
% rEnd      Coordinate of the end point.
% R         Radius of the arrow body.
% alpha     Angle of the head in degrees.
% lHead     Length of the head.
% N         Number of points on the curve, default value is stored in 
%           swpref.getpref('npatch').
%
% See also SWPLOT.CYLINDER.
%

if nargin == 0
    help swplot.arrow
    return
end

if nargin < 6
    N = swpref.getpref('npatch',[]);
end

rStart = rStart(:);
rEnd   = rEnd(:);
rArrow = rEnd-rStart;

alpha   = alpha * pi/180;
endHead = rStart+(1-lHead/norm(rEnd-rStart))*(rEnd-rStart);
rHead   = lHead * tan(alpha);

% generate two normal vectors perpendicular to the arrow direction
if any(cross(rArrow,[0; 0; 1]))
    nArrow1 = cross(rArrow,[0; 0; 1]);
else
    nArrow1 = cross(rArrow,[0; 1; 0]);
end
nArrow2 = cross(rArrow,nArrow1);
nArrow1 = nArrow1/norm(nArrow1);
nArrow2 = nArrow2/norm(nArrow2);

% generate points of a circle
phi     = linspace(0,2*pi,N);
cPoint  = R*(nArrow1*cos(phi)+nArrow2*sin(phi));

cPoint1 = bsxfun(@plus,cPoint,rStart);
cPoint2 = bsxfun(@plus,cPoint,endHead);
cPoint3 = bsxfun(@plus,cPoint*rHead/R,endHead);


V  = [cPoint1';cPoint2';cPoint3';rEnd(:)'];
L = (2:(N-2))';
% base circle
F1 = [ones(N-3,1) L L+1 nan(N-3,1)];
L  = (1:(N-1))';
% body cylinder
F2 = [L L+1 L+N+1 L+N];
% back of circle of the head
F3 = F2+N;
L  = (1:(N-1))'+2*N;
% head cone
F4 = [L L+1 ones(N-1,1)*(3*N+1) nan(N-1,1)];

F  = [F1;F2;F3;F4];
handle = patch('Vertices',V,'Faces',F,'FaceLighting','flat','EdgeColor','none','FaceColor','r');

end