function [M, k, n, name, pname, limit] = gm_planar(M0, x)
% planar magnetic structure constraint function 
%
% [M, k, n, name, pname, limit] = GM_PLANAR(M0, x) 
%
% It generates the parameters of arbitrary planar magnetic structure from
% phi angles, ordering wave vector and spin plane normal vector. All angles
% are in radian.
%
% Input:
%
% x         Input parameters in the following order:
%           (Phi1, Phi2, ... , kx, ky, kz, nTheta, nPhi).
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
% See also GM_SPHERICAL3D, GM_PLANARD.
%

if nargin == 0
    help gm_planar;
    return
end


if nargout <= 3
    
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
        M = bsxfunsym(@times,u*cos(phi) + v*sin(phi),M0);
    end
else
    nMagExt = size(M0,2);
    % provide the limits for the parameters
    name  = '2D planar structure';
    % parameter names
    pname = {};
    for ii = 1:nMagExt
        pname = [pname {sprintf('Phi%d_rad',ii)}]; %#ok<AGROW>
    end
    pname = [pname {'kx_rlu' 'ky_rlu' 'kz_rlu' 'nTheta_rad' 'nPhi_rad'}];
    % limits on input parameters
    limit = [zeros(1,nMagExt+5); [2*pi*ones(1, nMagExt) 1 1 1 pi 2*pi]];
    % garbage
    M = []; k = []; n = [];
end

end