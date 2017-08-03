function result = isop(symOp)
% function determines whether the matrix is a symmetry operator
%
% result = SWSYM.ISOP(Op)
%
% Input:
%
% Op        Symemtry operators with rotation and translation. Dimensions
%           are [3 4 nOp].
%
% See also SWSYM.GENERATOR, SWSYM.OPERATOR.
%

if nargin == 0
    help swsym.isop
    return
end

if size(symOp,1) == 3 && size(symOp,2) == 4
    result = true;
else
    result = false;
end

end