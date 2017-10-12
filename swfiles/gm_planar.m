function [S, k, n, name, pname, limit] = gm_planar(absS, x)
% planar magnetic structure constraint function 
% 
% ### Syntax
% 
% `[s, k, n, name, pname, limit] = gm_planar(S0, x)`
% 
% ### Description
% 
% `[s, k, n, name, pname, limit] = gm_planar(S0, x)` generates the
% parameters of arbitrary planar magnetic structure from $\varphi$ angles
% (in radian), ordering wave vector (rlu) and spin plane normal vector
% ($xyz$).
%  
% 
% ### Input Arguments
% 
% `x`
% : Input parameters in the following order: 
%   $[\varphi_1, \varphi_2, ... , k_x, k_y, k_z, n_\theta, n_\varphi]$.
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
% [gm_spherical3d] \| [gm_planard]
%
% *[rlu]: Reciprocal Lattice Unit
%

if nargin == 0
    help gm_planar;
    return
end


if nargout <= 3
    
    nMagExt = (length(x)-5);
    x       = x(:)';
    absS      = absS(:)';
    % Magnetic ordering wave vector in the crystallographic unit cell!
    k = x(end+(-4:-2));
    % Normal to the spin rotation plane.
    nTheta  = x(end-1);
    nPhi    = x(end);
    n = [sin(nTheta)*[cos(nPhi) sin(nPhi)] cos(nTheta)];
    % Angles in the spin plane.
    phi = x(1:nMagExt);
    
    [u, v] = sw_cartesian(n');
    
    if numel(absS)==1
        S = (u*cos(phi) + v*sin(phi))*absS;
    else
        % Check that the number of magnetic atoms is right
        if length(phi)~=length(absS)
            error('gm_planar:NumberOfMoments','The number of fitting parameters doesn''t produce the right number of moments!');
        end
        % Magnetic moments in orthogonal coordinate sysyem.
        S = bsxfun(@times,u*cos(phi) + v*sin(phi),absS);
    end
else
    nMagExt = size(absS,2);
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
    S = []; k = []; n = [];
end

end