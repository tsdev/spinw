function handle = cylinder(rStart, rEnd, R, N, close)
% draws a closed/open 3D cylinder
%
% [handle] = SWPLOT.CYLINDER(rStart, rEnd, R, {N}, {close})
%
% Input:
%
% rStart    Coordinate of the starting point.
% rEnd      Coordinate of the end point.
% R         Radius of the arrow body.
% N         Number of points on the curve, default value is stored in 
%           swpref.getpref('npatch').
% close     If true the cylinder is closed. Default is true.
%
% See also SWPLOT.ARROW.
%

if nargin == 0
    help swplot.cylinder
    return
end

if nargin < 4
    N = swpref.getpref('npatch',[]);
end

if nargin < 5
    close = true;
end

rArrow = rEnd(:) - rStart(:);
rStart = repmat(rStart(:),1,N);
rEnd   = repmat(rEnd(:),1,N);

% normal vectors to the cylinder axis
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
cPoint1 = cPoint+rStart;
cPoint2 = cPoint+rEnd;

% closing the begin of body
if close
    X = [rStart(1,:); cPoint1(1,:); cPoint2(1,:); rEnd(1,:)];
    Y = [rStart(2,:); cPoint1(2,:); cPoint2(2,:); rEnd(2,:)];
    Z = [rStart(3,:); cPoint1(3,:); cPoint2(3,:); rEnd(3,:)];
else
    X = [cPoint1(1,:); cPoint2(1,:)];
    Y = [cPoint1(2,:); cPoint2(2,:)];
    Z = [cPoint1(3,:); cPoint2(3,:)];   
end

handle = surface(X,Y,Z,'EdgeColor','none','FaceColor','r');

end