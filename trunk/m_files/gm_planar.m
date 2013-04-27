function [M, k, n] = gm_planar(M0, x)
% [M, k, n] = GM_PLANAR(M0, x) generates planar magnetic structure from phi
% angles, ordering wave vector and spin plane normal vector.
%
% x         Input parameters in the following order:
%           (Phi1, Phi2, ... , kx, ky, kz, nTheta, nPhi)
% M0        Size of magnetic moments: (M1, M2, ...) or scalar if all
%           moments are equal.
%
% See also GM_SPHERICAL3D.
%

nMagExt = (length(x)-5);
x       = x(:)';
M0      = M0(:)';
% Magnetic ordering wave vector in the crystallographic unit cell!
k = x(end+(-4:-2));
% Normal to the spin rotation plane.
nTheta  = x(end-1);
nPhi    = x(end);
n = [sin(nTheta)*[cos(nPhi) sin(nPhi)] cos(nTheta)];
% Angles in the spin plane.
phi = x(1:nMagExt);

[u, v] = sw_cartesian(n');

if numel(M0)==1
    M = (u*cos(phi) + v*sin(phi))*M0;
else
    % Check that the number of magnetic atoms is right
    if length(phi)~=length(M0)
        error('sw:gm_planar:NumberOfMoments','The number of fitting parameters doesn''t produce the right number of moments!');
    end
    % Magnetic moments in orthogonal coordinate sysyem.
    M = bsxfun(@times,u*cos(phi) + v*sin(phi),M0);
end
end