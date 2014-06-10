function type = sw_mattype(mat, epsilon)
% type = SW_MATTYPE(mat, {epsilon}) function determines the type of the mat
% 3x3xN matrix and returns a vector with dimensions of [1 N]. The following
% values are defined:
%
% Input:
%
% mat       Matrix with dimensions of [3 3 N].
% epsilon   Error bar on small matrix element. Defult is 1e-5.
%           Optional.
%
% Output:
%
% type      1   Heisenberg exchange
%           2   Anisotropic exchange
%           3   DM interaction
%           4   General matrix
%

if nargin==0
    help sw_mattype
    return
end

if nargin == 1
    epsilon = 1e-6*max(mat(:));
end

[mSize(1), mSize(2), mSize(3)] = size(mat);
type = zeros(1,mSize(3));

if any(mSize(1:2)-[3 3])
    error('sw:sw_mattype:InputError','Input matrix is not 3x3xN!');
end

if ~isreal(mat) && ~isa(mat,'sym')
    error('sw:sw_mattype:InputError','Input matrix is not real!');
end

for ii = 1:mSize(3)
    matI = mat(:,:,ii);
    dM = diag(matI);
    if sum(sum(abs(diag(dM)-matI))) < epsilon*6
        if sum(abs(dM-dM(1))) < epsilon*3
            typeT = 1;
        else
            typeT = 2;
        end
    elseif sum(sum(abs((matI+matI')))) < epsilon*6
        typeT = 3;
    else
        typeT = 4;
    end
    type(ii) = typeT;
end
end