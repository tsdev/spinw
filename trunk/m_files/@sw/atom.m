function atomList = atom(obj)
% generates all atomic positions in the unit cell
%
% atomList = ATOM(obj)
%
% Output:
% atomList is a structure with the following fields:
%
% r     Positions of the atoms in lattice units, dimensions are [3 nAtom].
% idx   Pointer to the atom in the unit_cell field, dimensions are
%       [nAtom 1].
% mag   Logical variable, whether the spin of the atom is non-zero,
%       dimensions are [nAtom 1].
%
% See also SW.MATOM.
%

% Generate all symmetry equivalent atoms
[atomList.r, atomList.idx] = sw_genatpos(obj.lattice.sym,obj.unit_cell.r);

atomList.mag = obj.unit_cell.S(atomList.idx)>0;

end % .atom