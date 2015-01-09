%% Bond
% Bonds are vectors connecting two magnetic ion.
%
%% Definition
% Bonds are defined on the crystal of magnetic ions using the positions
% within the unit cell. Each bond is defined by two atoms that are
% connected with direction. The bond points from atom1 to atom2 where atom1
% is at the zeroth unit cell (cell at the origin of the lattice coordinate
% system), atom2 is in the unit cell defined by the *dl* translation vector
% (in lattice units). The equivalent bonds that are related by lattice
% vectors to the above defined one are not stored explicitly but assumed in
% the spin wave calculation. The list of bonds are stored in the
% sw.coupling field, where each column of dl, atom1 and atom2 subfields
% defines different bonds. Each bond has an identifier stored in the idx
% subfield. These identifiers let the user select certain set of bonds
% easily. If two bonds have the same identifier, they are regarded
% equivalent. Coupling matrix can be quickly assigned to multiple
% equivalent bonds using the sw.addcoupling command.
%
%% Generating bonds
% Although the sw.coupling matrix can be filled manually, the
% sw.gencoupling command can generate the list of bonds automatically. The
% generated list of bonds sorted according to increasing length, however no
% particular order can be assumed between equal length bonds. If no
% symmetry operators are considered for the generation of bonds
% ('forceNoSym' option set to true) all bonds with equal length are
% assigned the same identifier, starting with 1 for the shortest bonds. If
% symmetry operators are considered, only symmetry equivalent bonds will
% have the same identifier. 
%
%% Listing bonds
% To list bonds in an easyer to read format, the sw.couplingtable command
% can be used. To list bonds with a set of identifiers use the
% sw.couplingtable(bond_id) command.






