function atomList = atom(obj)
% generates all atomic positions in the unit cell
%
% atomList = ATOM(obj)
%
% From the given atomic positions in the obj.unit_cell field the function
% generates all atomic positions using the symmetry operators of the given
% space group. If no symmetry is defined, the function returns simply the
% positions stored in obj.unit_cell.
%
% Input:
%
% obj   sw class object.
%
% Output:
%
% atomList is a structure with the following fields:
%
%   r       Positions of the atoms in lattice units, dimensions are 
%           [3 nAtom]. 
%   idx     Pointer to the atom in the unit_cell field, dimensions are
%           [nAtom 1].
%   mag     Logical variable, whether the spin of the atom is non-zero,
%           dimensions are [nAtom 1].
%
% Example:
%
% sw_addsym('x+1/2,y+1/2,z;x+1/2,y,z+1/2;x,y+1/2,z+1/2','FCC');
% 
% cryst = sw;
% cryst.genlattice('lat_const',[8 8 8],'sym','FCC')
% cryst.addatom('r',[0 0 0],'label','Atom1')
% atomList = cryst.atom;
%
% This will create a new symmetry, that contains all the translations of
% the FCC lattice. Then creates a crystal with an atom at [0 0 0] position.
% The cryst.atom lists all 4 symmetry equivalent positions generated using
% the 'FCC' symmetry operators.
%
% See also SW, SW.MATOM, SW_ADDSYM, SW.GENLATTICE, SW.ADDATOM.
%

% Generate all symmetry equivalent atoms
[atomList.r, atomList.idx] = sw_genatpos(obj.lattice.sym,obj.unit_cell.r);

if obj.symbolic
    atomList.mag = ~sw_always(obj.unit_cell.S(atomList.idx)==0);
else
    atomList.mag = obj.unit_cell.S(atomList.idx)>0;
end

end