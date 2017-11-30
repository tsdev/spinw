function [M, k, n, name, pname, limit] = gm_spherical3d(M0, x)
% magnetic structure constraint function with spherical parameterisation
% 
% ### Syntax
% 
% `[m, k, n, name, pname, limit] = gm_spherical3d(S0, x)`
% 
% ### Description
% 
% `[m, k, n, name, pname, limit] = gm_spherical3d(S0, x)` generates
% magnetic structure from given parameters while constraining the length of
% the spin on each atom. The parametrization of the magnetic structure
% consists of 2 spherical coordinates $(\theta,\varphi)$ angles per
% magnetic atom. All angles are in radian.
%   
% ### Input Arguments
%   
% `x`
% : Input parameters in the following order: 
%   $[\theta_1, \varphi_1, ... , k_x, k_y, k_z, n_\theta, n_\varphi]$.
% 
% `S0`
% : Spin quantum number in a row vector $(S_1, S_2, ...)$ or scalar if all
%   spins are equal.
%   
% ### Output Arguments
%   
% `S`
% : Matrix, containing the spin orientations with dimensions of $[3\times n_{magExt}]$.
%       Every column contains the $(S_x S_y S_z)$ spin components of
%       a magnetic atom in the $xyz$ coordinate system.
%
% `k`
% : Magnetic ordering wavevector in rlu units in a row vector.
%
% `n`
% : Normal vector around which the spins are rotating for non-zero
%       propagation vector in a row vector.
%
% `name`
% : String, storing the name of the function.
%
% `pname`
% : Name of the input parameters in a cell: `{'Phi1_rad', ...}`.
%
% `limit`
% : Limits on the input parameters in a matrix with dimensions of $[2\times n_{param}]$. Every
%       column contains a lower and upper limit on the corresponding
%       parameter.
%   
% ### See Also
%   
% [gm_planar]

if nargin == 0
    swhelp gm_spherical3d
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
            error('gm_spherical3d:NumberOfMoments','The number of fitting parameters doesn''t produce the right number of moments!');
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
