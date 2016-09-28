function [n, collinear] = sw_nvect(S, epsilon)
% determines the best normal vector for the set of vectors
%
% [n collinear] = SW_NVECT(S, {epsilon})
%
% The function can deal with complex vectors, separating the real and
% complex parts as separate vectors.
%
% S           Array of column vectors, dimensions are [3 N].
% epsilon     Upper limit of the collinearity and the lower limit of
%             coplanarity. Default value is 0.1. If epsilon = 1 the
%             function returns n vector closest to the collinear direction,
%             if epsilon = 2, n will be closest to the normal of the spin
%             plane.
%
% n           Vector parallel to the collinear spin direction or
%             perpendicular to the plane of the spins in planar structures,
%             dimensions are [1 3].
%
% collinear   If true, the set of vectors are collinear, then n is parallel to the
%             input vectors.
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