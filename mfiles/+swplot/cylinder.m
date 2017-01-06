function handle = cylinder(varargin)
% draws a closed/open 3D cylinder
%
% handle = SWPLOT.CYLINDER(rStart, rEnd, R, {N}, {close})
%
% handle = SWPLOT.CYLINDER(hAxis,...)
%
% draw to the specific axis
%
% Input:
%
% hAxis     Axis handle.
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

if numel(varargin{1}) == 1
    % first input figure handle
    hAxis = varargin{1};
    rStart  = varargin{2};
    rEnd    = varargin{3};
    R       = varargin{4};
    nArgExt = nargin-4;
    argExt  = {varargin{5:end}};
    
else
    hAxis = gca;
    rStart  = varargin{1};
    rEnd    = varargin{2};
    R       = varargin{3};
    nArgExt = nargin-3;
    argExt  = {varargin{4:end}};
end

if nArgExt > 0
    N = argExt{1};
else
    N = swpref.getpref('npatch',[]);
end

if nArgExt > 1
    close  = argExt{2};
else
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

handle = patch(hAxis,'Vertices',V,'Faces',F,'FaceLighting','flat',...
    'EdgeColor','none','FaceColor','r','Tag','cylinder');

end