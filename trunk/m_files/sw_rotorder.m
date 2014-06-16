function N = sw_rotorder(R, T)
% N = SW_ROTORDER(R, T) determines the order of the (R,T) symmetry
% operator, where R is a rotation matrix and T is a translation. Maximum
% order is 10 if the matrix is not a rotation matrix of any
% crystallographic point group.
%
% Input:
%
% R         Rotation matrix, dimensions are [3 3].
% T         Translation vector, dimensions are [3 1] optional.
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
    T = zeros(3,1);
end

tol = 1e-5;

N  = 1;
RN = R;
TN = T;

while ((norm(RN-eye(3))>tol) || norm(TN)>tol) && (N<10)
    RN = R*RN;
    TN = R*TN;
    N  = N + 1;
end

end