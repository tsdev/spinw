function A = mdiag(A,dim)
% Returns the diagonal elements along two selected dimensions.
%
% Input:
% A         Multidimensional input array.
% dim       Contains two number, that selects two dimensions.
%
% The diagonal is the diagonal of every A(... dim(1)...dim(2)...) matrix.
% The result is contracted along dim(2).
%
% The default value for dim is dim = [1 2].
%
% Examples:
% For 2D matrices mmat is identical to the Matlab built-in diag:
% A = [1 2; 3 4];
% C = mdiag(A)
%
% C will be [1;4].
%
% For multidimensional arrays:
% A = repmat([1 2; 3 4],[1 1 3]);
% B = mdiag(A,[1 2])
% B will be an array with dimensions of 2x3: [1 1 1;4 4 4].
%

if nargin == 0
    help mdiag;
    return;
end

if (nargin < 2)
    dim = [1 2];
end

if numel(dim)~=2
    error('mdiag:WrongInput','dim has to be a two element array!');
end


nD = ndims(A);

permIdx = [dim find(~ismember(1:nD,dim))];

A  = permute(A,permIdx);
sA = size(A);

B = repmat(eye(sA(1:2)),[1 1 sA(3:end)]);

A = sum(A.*B,2);

% permute back
permIdx2 = 1:nD;
permIdx2(permIdx) = permIdx2;
permIdx2 = [permIdx2 permIdx2(dim(2))];
permIdx2(dim(2)) = [];

A = permute(A,permIdx2);

end