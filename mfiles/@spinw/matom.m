function mAtomList = matom(obj)
% generates all magnetic atoms in the unit cell
%
% Same as spinw.atom, but only lists the magnetic atoms, which has non-zero
% spin.
%
% See also SPINW.ATOM.
%

if isempty(obj.matomstore)
    atomList      = obj.atom;
    
    mAtomList.r   = atomList.r(:,atomList.mag==1);
    mAtomList.idx = atomList.idx(:,atomList.mag==1);
    mAtomList.S   = obj.unit_cell.S(mAtomList.idx);
    
    obj.matomstore = mAtomList;
    
    % add listener to lattice and unit_cell fields
    obj.propl    = addlistener(obj,'lattice',  'PostSet',@obj.modmatom);
    obj.propl(2) = addlistener(obj,'unit_cell','PostSet',@obj.modmatom);
else
    mAtomList = obj.matomstore;
end

end