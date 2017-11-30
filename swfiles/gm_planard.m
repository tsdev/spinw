function [M, k, n, name, pname, limit] = gm_planard(M0, x)
% planar magnetic structure constraint function 
% 
% ### Syntax
% 
% `[m, k, n, name, pname, limit] = gm_planard(m0, x) `
% 
% ### Description
% 
% Same function as [gm_planar], except that the input angles are all in
% degree.
%  
%
% ### See Also
% 
% [gm_planar]
%

if nargin == 0
    swhelp gm_planard
    return
end


if nargout <= 3
    
    nMagExt = (length(x)-5);
    x       = x(:)';
    M0      = M0(:)';
    % Magnetic ordering wave vector in the crystallographic unit cell!
    k = x(end+(-4:-2));
    % Normal to the spin rotation plane.
    nTheta  = x(end-1)*pi/180;
    nPhi    = x(end)*pi/180;
    n = [sin(nTheta)*[cos(nPhi) sin(nPhi)] cos(nTheta)];
    % Angles in the spin plane.
    phi = x(1:nMagExt)*pi/180;
    
    [u, v] = sw_cartesian(n');
    
    if numel(M0)==1
        M = (u*cos(phi) + v*sin(phi))*M0;
    else
        % Check that the number of magnetic atoms is right
        if length(phi)~=length(M0)
            error('gm_planar:NumberOfMoments','The number of fitting parameters doesn''t produce the right number of moments!');
        end
        % Magnetic moments in orthogonal coordinate sysyem.
        M = bsxfun(@times,u*cos(phi) + v*sin(phi),M0);
    end
else
    nMagExt = size(M0,2);
    % provide the limits for the parameters
    name  = '2D planar structure';
    % parameter names
    pname = {};
    for ii = 1:nMagExt
        pname = [pname {sprintf('Phi%d_deg',ii)}]; %#ok<AGROW>
    end
    pname = [pname {'kx_rlu' 'ky_rlu' 'kz_rlu' 'nTheta_deg' 'nPhi_deg'}];
    % limits on input parameters
    limit = [zeros(1,nMagExt+5); [360*ones(1, nMagExt) 1 1 1 180 360]];
    % garbage
    M = []; k = []; n = [];
end

end