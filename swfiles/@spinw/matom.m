function mAtomList = matom(obj)
% generates all magnetic atoms in the unit cell
%
% mAtomList = MATOM(obj)
%
% Same as spinw.atom, but only lists the magnetic atoms, which has non-zero
% spin. Also this function stores the generated list in spinw.cache.
%
% Output:
%
% mAtomList is a structure with the following fields:
%   r       Position of the magnetic atoms in a matrix with dimensions of 
%           [3 nMagAtom].
%   idx     Index in the symmetry inequivalent atom list spinw.unit_cell 
%           stored in a row vector with dimensions of [1 nMagAtom].
%   S       Spin of the magnetic atoms stored in a row vector with 
%           dimensions of [1 nMagAtom].
%
% See also SPINW.ATOM.
%

if isempty(obj.cache.matom)
    atomList      = obj.atom;
    
    mAtomList.r   = atomList.r(:,atomList.mag==1);
    mAtomList.idx = atomList.idx(:,atomList.mag==1);
    mAtomList.S   = obj.unit_cell.S(mAtomList.idx);
    
    obj.cache.matom = mAtomList;
    
    % add listener to lattice and unit_cell fields
    obj.addlistenermulti(1);
else
    mAtomList = obj.cache.matom;
end

end