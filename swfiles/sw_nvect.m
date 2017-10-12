function [n, collinear] = sw_nvect(S, epsilon)
% determines the best normal vector for the set of vectors
% 
% ### Syntax
% 
% `[n, collinear] = sw_nvect(V)`
% 
% `[n, collinear] = sw_nvect(V,epsilon)`
%
% ### Description
% 
% `[n, collinear] = sw_nvect(V)` determines whether the given set of
% vectors are collinear or coplanar. If they are coplanar, it returns the
% best fitting normal vector, while if they are collinear returns the
% average of the given vector.
%
% The function can also deal with complex vectors, separating the real and
% complex parts as separate vectors.
%
%
% `[n, collinear] = sw_nvect(V,epsilon)` also gives the upper limit of the
% collinearity controlled by `epsilon`.
%  
% ### Input Arguments
%
% `V`
% : Matrix of column vectors with dimensions of $[3\times N]$. Where each
%   column defines a vector.
%
% `epsilon`
% : Defines the limits of collinearity with the following values:
%   * `1`   the function always return the `n` closest to the collinear
%           direction,
%   * `2`   the function always return the `n` vector closest to the normal 
%           of the coplanar plane.
%   * `e`   upper limit of collinearity, default value is 0.1, smaller
%           positive values mean stricter limits on collinearity.
%  
% ### Output Arguments
%
% `n`
% : Row vector parallel to the collinear vector direction or
%   perpendicular to the best fitting plane of the coplanar vectors.
%  
% `collinear`
% : `true` if the given set of vectors are collinear.
%  

if nargin == 0
    help sw_nvect
    return
end

if nargin == 1
    epsilon = 0.1;
end

if ~isreal(S)
    % deal with complex vectors
    S = [real(S) imag(S)];
    S(:,sum(S.^2,1)==0) = [];
end

Srot   = bsxfun(@rdivide,S,sqrt(sum(S.^2,1)));
SS     = sum(repmat(permute(Srot,[1 3 2]),[1 3 1]).*repmat(permute(Srot,[3 1 2]),[3 1 1]),3);
[V, D] = eig(SS/size(Srot,2));
D      = diag(D);

% Fix result
if epsilon == 1
    collinear = true;
    n = V(:,3)';
elseif epsilon == 2
    collinear = false;
    n = V(:,1)';
    % D contains the eigenvectors of SS in increasing order.
    % SS is positive definit, thus all eigenvectors are positive.
elseif D(2) < epsilon
    % The structure is rather collinear.
    % The starting vector, size (1,3):
    collinear = true;
    n = V(:,3)';
else
    % The structure is rather planar or 3D.
    collinear = false;
    n = V(:,1)';
end

end