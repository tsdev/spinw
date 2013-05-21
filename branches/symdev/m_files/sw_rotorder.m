function N = sw_rotorder(R, tol)
% N = SW_ROTORDER(R, {tol}) determines the order of the symOp rotation matrix.
% Maximum is 10 if the matrix is not a rotation matrix of any
% crystallographic point group.
%
% Input:
%
% R         Rotation matrix, dimensions are [3 3].
% tol       Tolerance, default is 1e-5.
%
% Example:
% R^sw_order(R) == eye(3);
%
% See also SW_GENATPOS, SW_BASISMAT.
%

if nargin == 0
    help sw_rotorder;
    return;
end

if nargin == 1
    tol = 1e-5;
end

N = 1;

while (norm(R^N-eye(3))>tol) && (N<10)
    N = N + 1;
end

end