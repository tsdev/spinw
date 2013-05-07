function pOp = sw_pointsym(sym, r)
% pOp = SW_POINTSYM(sym, r) determines point group symmetry in an arbitrary
% position in the unit cell in any space group.
%
% Input:
%
% sym           Line index in the symmetry.dat file or string of the
%               symmetry operators or a cell containing the output of
%               sw_gencoord().
% r             Position in the unit cell, dimensions are [3 1].
%
% Output:
%
% pOp           Point group operators, dimensions are [3 3 npOp], these
%               operators act on the relative atomic positions (they are in
%               the lattice coordinate system). To convert them to
%               Cartesian coordinate system, use:
%                   R = A*pOp(:,:,ii)*inv(A)
%               Where A is a 3x3 matrix, containing the basis vectors of
%               the lattice as column vectors.
%
% See also SW_GENATPOS, SW_GENCOORD, SW_GENSYM, SW_BASISMAT.
%

if nargin == 0
    help sw_pointsym;
    return;
end

if ~iscell(sym)
    [symOp, symTr] = sw_gencoord(sym);
else
    symOp = sym{1};
    symTr = sym{2};
end

[~, ~, isMoved] = sw_genatpos({symOp symTr}, r);

% point group operators are the ones that does NOT move the atom
pOp = symOp(:,:,~isMoved{1});

end