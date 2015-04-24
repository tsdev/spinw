function hOrbital = sw_t2g(idx, varargin)
% plots t2g electron orbital surfaces
%
% hOrbital = SW_T2G(idx, 'Option1',Value1,...)
%
% Input:
%
% idx   Index of the orbitaL:
%           1 d_xy
%           2 d_xz
%           3 d_yz
%
% Options:
%
% r0            Origin coordinates of the orbital, size is (1,3), default
%               is [0 0 0].
% v1            Vector pointing towards the first positive charge, defines
%               the x-axis, size is (1,3), default is [1 0 0].
% v2            Vector pointing towards the second positive charge, defines
%               the y-axis, size is (1,3), default is [0 1 0].
% surfRes       Number of points in the surface mesh along every axis,
%               default is 30.
% rLim          Limits of the axes, default is 32.
% P             Constant probability surface, default is 1E-5.
% rBohr         Bohr radius, default is 1.
% norm          Whether to normalise the axes, default is true.
% scale         Scale of the radius of the orbital, default is 1.
% plotv         Whether to plot the v1, v2 vectors. Default is false.
%
% See also SW_ORBITAL, SW_DRAWPOLY, SW_ADDOBJECT.
%

if nargin == 0
    help sw_t2g;
    return;
end

inpForm.fname  = {'r0'    'v1'    'v2'    'scale' 'plotv'};
inpForm.defval = {[0 0 0] [1 0 0] [0 1 0] 1       false  };
inpForm.size   = {[1 3]   [1 3]   [1 3]   [1 1]   [1 1]  };

param = sw_readparam(inpForm,varargin{:});

v1 = param.v1;
v2 = param.v2;

n = 3;
l = 2;

switch idx
    case 1
        m = 2;
        osign = 2;
    case 2
        m = 1;
        osign = 1;
    case 3
        m = 1;
        osign = 2;
    otherwise
        error('sw:sw_t2g:WrongIdx','Idx has to be 1,2 or 3!');
end

% create surface polyeder of the orbital
surfO = sw_orbital([n l m osign], varargin{:});

% rotate
% base plane made by two vectors

% rotation of xy plane to the v1-v2 plane
n = cross(v1,v2);
rotAxis  = cross([0 0 1],n);
rotAngle = acos([0 0 1]*n'/norm(n));

if abs(rotAngle) > 1e-5
    [~, R] = sw_rot(rotAxis,rotAngle);
    surfO.vertices = surfO.vertices*R';
else
    R = eye(3);
end

% rotation of the x' axis to v1
xp = [1 0 0]*R';
rotAxis2  = cross(xp,v1);
rotAngle2 = acos(xp*v1'/norm(v1));

if abs(rotAngle2) > 1e-5
    % rotate away from the ions
    [~, R2] = sw_rot(rotAxis2,rotAngle2);
    surfO.vertices = surfO.vertices*R2';
end

% scale the orbital
surfO.vertices = surfO.vertices * param.scale;

% shift the origin to r0
surfO.vertices = bsxfun(@plus,surfO.vertices,param.r0(:)');

hOrbital = patch(surfO);

set(hOrbital,'FaceColor',[1 0 0],'FaceAlpha',1,'FaceLighting','phong','EdgeColor', 'none'); %'FaceColor','interp'

if param.plotv
    r0 = param.r0;
    l1 = line(r0(1)+[0 v1(1)],r0(2)+[0 v1(2)],r0(3)+[0 v1(3)],'linewidth',5,'color',[0 0 0]);
    l2 = line(r0(1)+[0 v2(1)],r0(2)+[0 v2(2)],r0(3)+[0 v2(3)],'linewidth',5,'color',[0 0 0]);
    
    hOrbital(2:3) = [l1 l2];
end

end