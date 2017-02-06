function [M, k, n, name, pname, limit] = gm_spherical3d(M0, x)
% magnetic structure constraint function with spherical parameterisation
%
% [M, k, n, name, pname, limit] = GM_SPHERICAL3D(M0, x) 
%
% It generates the parameters of magnetic moments and normal vector from
% spherical (theta,phi) coordinates. All angles are in radian.
%
% Input:
%
% x         Input parameters in the following order:
%           (Theta1, Phi1, Theta2, Phi2, ... , kx, ky, kz, nTheta, nPhi).
% M0        Size of magnetic moments: (M1, M2, ...) or scalar if all
%           moments are equal.
%
% Output:
%
% M         Array, containing the magnetic moments, dimensions are
%           [3 nMagExt]. Every column contain the [Mx; My; Mz] magnetic
%           moment components of a magnetic atom in the xyz coordinate
%           system.
% k         Magnetic ordering wavevector in r.l.u., dimensions are [1 3].
% n         Normal vector to the plane of the incommensurate spins (if k
%           non-zero).
%
% Optional outputs:
% only produced if the output is requested.
%
% name      Name of the function.
% pname     Name of the input parameters in a cell: {'Param1' 'Param2',...}
% limit     Default limits on the input parameters, dimensions are [2 nX].
%           Every column contains a lower and upper limit on the parameter.
%
% See also GM_PLANAR.
%

if nargin == 0
    help gm_spherical3d;
    return
end

if nargout <= 3
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
    end
    
    % Normal to the spin rotation plane.
    nTheta  = x(end-1);
    nPhi    = x(end);
    n = [sin(nTheta)*[cos(nPhi) sin(nPhi)] cos(nTheta)];
    % Magnetic ordering wave vector in the crystallographic unit cell!
    k = x(end+(-4:-2));
    
else
    nMagExt = size(M0,2);
    % provide the limits for the parameters
    name  = '3D structure with spherical coordinates';
    % parameter names
    pname = {};
    for ii = 1:nMagExt
        pname = [pname {sprintf('Theta%d_rad',ii) sprintf('Phi%d_rad',ii)}]; %#ok<AGROW>
    end
    pname = [pname {'kx_rlu' 'ky_rlu' 'kz_rlu' 'nTheta_rad' 'nPhi_rad'}];
    % limits on input parameters
    limit = [zeros(1,nMagExt*2+5); [repmat([pi 2*pi],[1 nMagExt]) 1 1 1 pi 2*pi]];
    % garbage
    M = []; k = []; n = [];
end

end
