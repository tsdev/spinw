function [handle] = circle(varargin)
% creates a circle surface in 3 dimensions
%
% handle = SWPLOIT.CIRCLE(r0, n, R, {N})
%
% handle = SWPLOIT.CIRCLE(hAxis,...)
%
% plot to specific axis.
%
% Input:
%
% hAXis     Axis handle.
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

if numel(varargin{1}) == 1
    % first input figure handle
    hAxis   = varargin{1};
    r0      = varargin{2};
    n       = varargin{3};
    R       = varargin{4};
    nArgExt = nargin-4;
    argExt  = {varargin{5:end}};
    
else
    hAxis   = gca;
    r0      = varargin{1};
    n       = varargin{2};
    R       = varargin{3};
    nArgExt = nargin-3;
    argExt  = {varargin{4:end}};
end

if nArgExt > 0
    N = argExt{1};
else
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

handle = patch(hAxis,edge{:},'FaceLighting','flat','EdgeColor','none',...
    'FaceColor','r','Tag','circle');

end