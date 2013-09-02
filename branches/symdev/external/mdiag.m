function C = mdiag(A,dim)
% Returns the diagonal elements along two selected dimensions.
%
% Input:
% A         Multidimensional input array.
% dim       Contains two number, that selects two dimensions.
%
% The diagonal is the diagonal of every A(... dim(1)...dim(2)...) matrix.
% The dim(2) dimension will be contracted.
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

sA = size(A);
nD = ndims(A);


sM = sA(dim);
sR = sA; sR(dim) = 1;

pVect = 1:nD; pVect = pVect()
repmat(permute(eye(sM),pVect())







nA = [size(A),ones(1,nD-nDA)]; nA = nA(dim); 

% form A matrix
% (nA1) x (nA2) x nB2
A = repmat(A,[ones(1,nD) nB(2)]);
% form B matrix
% nA1 x (nB1) x (nB2)
idx = 1:nD+1; idx([dim end]) = idx([end dim]);
repv = ones(1,nD+1); repv(dim(1)) = nA(1);

B = repmat(permute(B,idx),repv);

% multiply with expanding along singleton dimensions
C = sum(bsxfun(@times,A,B),dim(2));

idx2 = 1:nD+1; idx2([dim end]) = idx2([dim(1) end dim(2)]);

% permute back the final result to the right size
C = permute(C,idx2);

end