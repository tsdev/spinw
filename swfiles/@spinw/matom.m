function mAtomList = matom(obj)
% generates magnetic lattice
% 
% ### Syntax
% 
% `mAtomList = matom(obj)`
% 
% ### Description
% 
% `mAtomList = matom(obj)` is the same as [spinw.atom], but only lists the
% magnetic atoms, which have non-zero spin. Also this function stores the
% generated list in [spinw.cache].
% 
% ### Output Arguments
% 
% `mAtomList`
% : structure with the following fields:
%   * `r`   Position of the magnetic atoms in a matrix with dimensions of 
%     $[3\times n_{magAtom}]$.
%   * `idx` Index in the symmetry inequivalent atom list [spinw.unit_cell] 
%     stored in a row vector with $n_{magAtom}]$ number of elements.
%   * `S`   Spin of the magnetic atoms stored in a row vectorwith
%     $n_{magAtom}]$ number of elements.
% 
% ### See Also
% 
% [spinw.atom]
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