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

% Defines the number of independent atoms.
nAtom0 = size(obj.unit_cell.r,2);

% It stores parameters of all atoms.
atomList.r   = zeros(3,0); % (3,nAtom)
atomList.idx = zeros(0,1); % (nAtom,1)
atomList.mag = zeros(0,1); % (nAtom,1)

% It stores atomic paramters in atom variable.
for ii = 1:nAtom0
    
    rNew = sw_genatpos(obj.lattice.sym, obj.unit_cell.r(:,ii));
    
    if obj.unit_cell.S(ii)>0
        atomList.mag = [atomList.mag ones(1,size(rNew,2)) ];
    else
        atomList.mag = [atomList.mag zeros(1,size(rNew,2))];
    end
    
    atomList.r   = [atomList.r rNew];
    atomList.idx = [atomList.idx repmat(ii,1,size(rNew,2))];
end

end % .atom
