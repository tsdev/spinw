function [M, k, n] = gm_spherical3d(M0, x)
% [M, k, n] = GM_SPHERICAL3D(M0, x) generates magnetic moments and normal
% vector from spherical (theta,phi) coordinates.
%
% x         Input parameters in the following order:
%           (Theta1, Phi1, Theta2, Phi2, ... , kx, ky, kz, nTheta, nPhi)
% M0        Size of magnetic moments: (M1, M2, ...) or scalar if all
%           moments are equal.
%
% See also GM_PLANAR.
%

nMagExt = (length(x)-5)/2;
x       = x(:)';
M0      = M0(:)';
MTheta  = x((1:nMagExt)*2-1);
MPhi    = x((1:nMagExt)*2);

if numel(M0) == 1
    % Magnetic moments in orthogonal coordinate sysyem.
    M = [sin(MTheta).*cos(MPhi); sin(MTheta).*sin(MPhi); cos(MTheta)]*M0;
else
    % Check that the number of magnetic atoms is right
    if length(MTheta)~=length(M0)
        error('sw:gm_spherical3d:NumberOfMoments','The number of fitting parameters doesn''t produce the right number of moments!');
    end
    % Magnetic moments in orthogonal coordinate sysyem.
    M = bsxfun(@times,[sin(MTheta).*cos(MPhi); sin(MTheta).*sin(MPhi); cos(MTheta)],M0);
    % Magnetic ordering wave vector in the crystallographic unit cell!
    k = x(end+(-4:-2));
    % Normal to the spin rotation plane.
    n = x(end+(-1:0));
end
end
