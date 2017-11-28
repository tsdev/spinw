function N = oporder(symOp)
% determine the order of the symmetry operator
% 
% ### Syntax
% 
% `N = swsym.oporder(symOp)`
% 
% ### Description
% 
% `N = swsym.oporder(symOp)` determines the order of the `symOp` symmetry
% operator, where `symOp(:,1:3)` is a rotation matrix and `symOp(:,4)` is a
% translation. The value of 10 is returned if the matrix is not a valid
% crystallographic symmetry operator.
% 
% ### Examples
% 
% Raising any operator to the calculated order will alway return identity:
%
% ```
% >>O = swsym.generator('y,z,x')>>
% >>R = O(:,1:3)^swsym.oporder(O)>>
% ```
% 
% ### Input Arguments
% 
% `symOp`
% :	Symmetry operator in a matrix with dimensions of $[3\times 4]$.
% 
% ### Output Arguments
%
% `N`
% : Integer, the order of the operator.
%
% ### See Also
% 
% [swsym.generator] \| [sw_basismat]
%

if nargin == 0
    help swsym.oporder
    return
end

N  = 1;
RN = symOp(:,1:3);
TN = round(symOp(:,4)*12);

while (norm(RN-eye(3))>10*eps || norm(TN)) && (N<10)
    RN = symOp(:,1:3)*RN;
    TN = mod(round(symOp(:,1:3)*TN+symOp(:,4)),12);
    N  = N + 1;
end

end