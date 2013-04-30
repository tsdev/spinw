function type = sw_mattype(mat)
% type = SW_MATTYPE(mat) function determines the type of the mat 3x3xN
% matrix and returns a vector with dimensions of [1 N]. The following
% values are defined:
%
% type    1   Heisenberg exchange
%         2   Anisotropic exchange
%         3   DM interaction
%         4   General matrix
%

if nargin==0
    help sw_mattype
    return
end

[mSize(1), mSize(2), mSize(3)] = size(mat);
type = zeros(1,mSize(3));

if any(mSize(1:2)-[3 3])
    error('sw:sw_mattype:InputError','Input matrix is not 3x3xN!');
end

if ~isreal(mat)
    error('sw:sw_mattype:InputError','Input matrix is not real!');
end

for ii = 1:mSize(3)
    matI = mat(:,:,ii);
    dM = diag(matI);
    if ~any(any(diag(dM)-matI))
        if ~any(dM-dM(1))
            typeT = 1;
        else
            typeT = 2;
        end
    elseif ~any(any(matI+matI'))
        typeT = 3;
    else
        typeT = 4;
    end
    type(ii) = typeT;
end
end