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

V = [cPoint1';cPoint2'];

% generate faces
L = (1:(N-1))';
F = [L L+1 L+N+1 L+N];

if close
    % caps
    L = (2:(N-2))';
    F1 = [ones(N-3,1) L L+1 nan(N-3,1)];
    F2 = F1+N;
    F = [F1;F;F2];
end

handle = patch('Vertices',V,'Faces',F,'FaceLighting','flat','EdgeColor','none','FaceColor','r');

end