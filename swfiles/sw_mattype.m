function type = sw_mattype(mat, epsilon)
% classifies square matrices
% 
% ### Syntax
% 
% `type = sw_mattype(mat)`
% 
% `type = sw_mattype(mat,epsilon)`
%
% ### Description
% 
% `type = sw_mattype(mat)` determines the type of the input matrix `mat`
% which stacked $[3\times 3]$ matrices. It determines the type of exchnge
% interaction that the matrix belongs to.
% 
% {{note Also works on symbolic matrices, but keep all symbols real for consistent
% result!}}
%
% ### Input Arguments
% 
% `mat`
% : Matrix with dimensions of $[3\times 3\times N]$.
% 
% `epsilon`
% : optional error bar on small matrix elements, default value is $10^{-5}$.
% 
% ### Output Arguments
% 
% `type`
% : Row vector with $N$ elements each having one of the following value:
%   * `1`   Heisenberg exchange,
%   * `2`   anisotropic exchange,
%   * `3`   DM interaction,
%   * `4`   general matrix.
%
% *[DM]: Dzyaloshinskii-Moriya
%

if nargin==0
    swhelp sw_mattype
    return
end

if nargin == 1 && ~isa(mat,'sym')
    epsilon = max(1e-6*max(mat(:)),1e-10);
end

[mSize(1), mSize(2), mSize(3)] = size(mat);
type = zeros(1,mSize(3));

if any(mSize(1:2)-[3 3])
    error('sw_mattype:InputError','Dimensions of the Input matrix is not 3x3xN!');
end

if ~isreal(mat) && ~isa(mat,'sym')
    error('sw_mattype:InputError','Input matrix is not real or symbolic!');
end

if ~isa(mat,'sym')
    % double type matrix
    for ii = 1:mSize(3)
        matI = mat(:,:,ii);
        dM = diag(matI);
        if sum(sum(abs(diag(dM)-matI))) <= epsilon*6
            if sum(abs(dM-dM(1))) <= epsilon*3
                % unity matrix * scalar
                typeT = 1;
            else
                % anisotropic diagonal
                typeT = 2;
            end
        elseif sum(sum(abs((matI+matI')))) <= epsilon*6
            % pure antisymmetric
            typeT = 3;
        else
            % general matrix
            typeT = 4;
        end
        type(ii) = typeT;
    end
else
    % symbolic matrix
    % keep symbolic variables all real for consistent results
    for ii = 1:mSize(3)
        matI = mat(:,:,ii);
        dM = diag(matI);
        if logical(sum(sum(abs(diag(dM)-matI))) == 0)
            if logical(sum(abs(dM-dM(1))) == 0)
                % unity matrix * scalar
                typeT = 1;
            else
                % anisotropic diagonal
                typeT = 2;
            end
        elseif logical(sum(sum(abs((matI+matI')))) == 0)
            % pure antisymmetric
            typeT = 3;
        else
            % general matrix
            typeT = 4;
        end
        type(ii) = typeT;
    end
end

end