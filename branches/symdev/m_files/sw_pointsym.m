function pOp = sw_pointsym(sym, r)
% pOp = SW_POINTSYM(sym, r) determines point group symmetry in an arbitrary
% position in the unit cell in any space group.
%
% Input:
%
% sym           Line index in the symmetry.dat file or string of the
%               symmetry operators.
% r             Position in the unit cell, dimensions are [3 1].
%
% Output:
%
% pOp           Point group operators, dimensions are [3 3 npOp].
%
% See also SW_GENATPOS, SW_GENCOORD, SW_GENSYM, SW_BASISMAT.
%

if nargin == 0
    help sw_pointsym;
    return;
end

[symOp, ~] = sw_gencoord(sym);

[~, ~, isMoved] = sw_genatpos(sym, r);

% point group operators are the ones that does NOT move the atom
pOp = symOp(:,:,~isMoved{1});

end