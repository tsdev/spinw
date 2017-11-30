function [M, k, n, name, pname, limit] = gm_spherical3dd(M0, x)
% magnetic structure constraint function with spherical parameterisation
% 
% ### Syntax
% 
% `[m, k, n, name, pname, limit] = gm_spherical3dd(m0, x) `
% 
% ### Description
% 
% Same function as [gm_spherical3d], except that the input angles are all in
% degree.
%  
%
% ### See Also
% 
% [gm_spherical3d]
%


if nargin == 0
    swhelp gm_spherical3dd
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
        M = [sind(MTheta).*cosd(MPhi); sind(MTheta).*sind(MPhi); cosd(MTheta)]*M0;
    else
        % Check that the number of magnetic atoms is right
        if length(MTheta)~=length(M0)
            error('gm_spherical3d:NumberOfMoments','The number of fitting parameters doesn''t produce the right number of moments!');
        end
        % Magnetic moments in orthogonal coordinate sysyem.
        M = bsxfun(@times,[sind(MTheta).*cosd(MPhi); sind(MTheta).*sind(MPhi); cosd(MTheta)],M0);
    end
    
    % Normal to the spin rotation plane.
    nTheta  = x(end-1);
    nPhi    = x(end);
    n = [sind(nTheta)*[cosd(nPhi) sind(nPhi)] cosd(nTheta)];
    % Magnetic ordering wave vector in the crystallographic unit cell!
    k = x(end+(-4:-2));
    
else
    nMagExt = size(M0,2);
    % provide the limits for the parameters
    name  = '3D structure with spherical coordinates';
    % parameter names
    pname = {};
    for ii = 1:nMagExt
        pname = [pname {sprintf('Theta%d_deg',ii) sprintf('Phi%d_deg',ii)}]; %#ok<AGROW>
    end
    pname = [pname {'kx_rlu' 'ky_rlu' 'kz_rlu' 'nTheta_deg' 'nPhi_deg'}];
    % limits on input parameters
    limit = [zeros(1,nMagExt*2+5); [repmat([180 360],[1 nMagExt]) 1 1 1 180 360]];
    % garbage
    M = []; k = []; n = [];
end

end
