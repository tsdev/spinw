function handle = arrow(rStart, rEnd, R, alpha, lHead, N)
% draws a 3D arrow
%
% [handle] = SWPLOT.ARROW(rStart, rEnd, R, alpha, lHead, N)
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
    help swplot.arrow
    return
end

nan2 = nan(2,1);

rStart = rStart(:);
rEnd   = rEnd(:);
rArrow = rEnd-rStart;

alpha   = alpha * pi/180;
endHead = rStart+(1-lHead/norm(rEnd-rStart))*(rEnd-rStart);
rHead   = lHead * tan(alpha);

%body
if any(cross(rArrow,[0; 0; 1]))
    nArrow1 = cross(rArrow,[0; 0; 1]);
else
    nArrow1 = cross(rArrow,[0; 1; 0]);
end

nArrow2 = cross(rArrow,nArrow1);
nArrow1 = nArrow1/norm(nArrow1);
nArrow2 = nArrow2/norm(nArrow2);
phi     = linspace(0,2*pi,N);
cPoint  = R*(nArrow1*cos(phi)+nArrow2*sin(phi));
cPoint1 = bsxfun(@plus,cPoint,rStart);
cPoint2 = bsxfun(@plus,cPoint,endHead);

% closing the begin of body
cPoint3 = bsxfun(@plus,cPoint,rStart);
X = [cPoint3(1,:); rStart(1)*ones(1,N)];
Y = [cPoint3(2,:); rStart(2)*ones(1,N)];
Z = [cPoint3(3,:); rStart(3)*ones(1,N)];

X = [X nan2 [cPoint1(1,:); cPoint2(1,:)]];
Y = [Y nan2 [cPoint1(2,:); cPoint2(2,:)]];
Z = [Z nan2 [cPoint1(3,:); cPoint2(3,:)]];

% annulus for head
cPoint4 = bsxfun(@plus,cPoint*rHead/R,endHead);
X = [X nan2 [cPoint4(1,:); endHead(1)*ones(1,N)]];
Y = [Y nan2 [cPoint4(2,:); endHead(2)*ones(1,N)]];
Z = [Z nan2 [cPoint4(3,:); endHead(3)*ones(1,N)]];

% cone of head
X = [X nan2 [cPoint4(1,:); rEnd(1)*ones(1,N)]];
Y = [Y nan2 [cPoint4(2,:); rEnd(2)*ones(1,N)]];
Z = [Z nan2 [cPoint4(3,:); rEnd(3)*ones(1,N)]];

handle = surface(X,Y,Z,'LineStyle','none','EdgeColor','none','FaceColor','r');

end