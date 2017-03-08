function N = oporder(symOp)
% determine the order of the symmetry operator
%
% N = SWSYM.OPORDER(symOp)
%
% It determines the order of the symOp symmetry operator, where
% symOp(:,1:3) is a rotation matrix and symOp(:,4) is a translation.
% Maximum order is 10 if the matrix is not a rotation matrix of any
% crystallographic point group.
%
% Input:
%
% symOp 	Symmetry operator in a matrix.
%
% Example:
%
% R^sw_symorder([R zeros(3,1)]) == eye(3);
%
% See also SWSYM.GENERATOR, SW_BASISMAT.
%

if nargin == 0
    help swsym.oporder
    return
end

N  = 1;
RN = symOp(:,1:3);
switch size(symOp,2)
    case 3
        symOp(:,4) = zeros(3,1);
    case 4
        % do nothing
    otherwise
        error('oporder:WrongInput','The given matrix is not a symmetry operator!')
end

TN = round(symOp(:,4)*12);

while (norm(RN-eye(3))>10*eps || norm(TN)) && (N<10)
    RN = symOp(:,1:3)*RN;
    TN = mod(round(symOp(:,1:3)*TN+symOp(:,4)),12);
    N  = N + 1;
end

end