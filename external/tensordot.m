function A = tensordot(A,B,dim1,dim2)
% Tensor product of arbitrary dimensional matrices.
%
% C = tensorprod(A,B,dimA,{dimB})
%
% Input:
% A, B      Multidimensional input arrays.
% dim       Contains two number, that selects two dimensions.
%
% The multiplication is a standard matrix multiplication. The two matrices
% are selected by dim:
%   AB = A(... dim(1) ... dim(2) ...) * B(... dim(1) ... dim(2) ...)
% The necessary condition that the multiplication can be performed:
%   size(A,dim(2)) = size(B,dim(1))
%
% Singleton dimensions in both A and B matrices are supported.
%
% The default value for dim is dim = [1 2].
%
% Examples:
% For 2D matrices mmat is identical to the Matlab built-in multiplication:
% A = [1 2 3 4];
% B = [1;2;3;4];
% C = mmat(A,B)
%
% C will be 30.
%
% For multidimensional arrays:
% A = repmat([1 2 3 4],[1 1 5]);
% B = [1 2 3 4]';
% C = mmat(A,B)
% C will be an array with dimensions of 1x1x5 and every element is 30.
%

if nargin == 0
    help tensordot;
    return;
end

if nargin < 3
    dim1 = [1 2];
end
if nargin < 4
    dim2 = dim1;
end

if numel(dim1)~=numel(dim2)
    error('tensordot:WrongInput',['The number of indices to sum has to be '...
        'the same on both matrix!']);
end
% check the sizes of matrices
for ii = 1:numel(dim1)
    if size(A,dim1(ii))~=size(B,dim2(ii)) && (size(A,dim1(ii))~=1 && size(B,dim2(ii))~=1)
        error('tensordot:WrongInput',['Input matrices have '...
            'incompatible length along the %dth dimension!'],ii)
    end
end

% permute the matrices find dimensions to put to the left
nA = ndims(A);
nB = ndims(B);

sA0 = size(A); 
sB0 = size(B);
% add extra one up to max(dim1)
if max(dim1)>nA
    sA0 = [sA0 ones(1,max(dim1)-nA)];
    nA = max(dim1);
end
if max(dim2)>nB
    sB0 = [sB0 ones(1,max(dim2)-nB)];
    nB = max(dim2);
end

dA = 1:nA; 
dB = 1:nB;

dA(dim1) = [];
dB(dim2) = [];

sD1 = sA0(dim1);
sD2 = sB0(dim2);

sA = sA0; sA(dim1) = [];
sB = sB0; sB(dim2) = [];

% for example A(i,j,k,l) --> A((jk),(il)) for dim1 = [1 4]
A = permute(A,[dA dim1]);
B = permute(B,[dim2 dB]);
% reshape the matrix
A = reshape(A,[prod(sA) prod(sD1)]);
B = reshape(B,[prod(sD2) prod(sB)]);

% 2 step procedure
% #1 multiply non-singleton dimensions

% #2 multiply singleton dimensions


if any(sD1 == 1)
    % add expanded dimensions
    A = reshape(bsxfun(@times,A,B.'),[sA sB]);
else
    A = reshape(A*B,[sA sB]);
end

end