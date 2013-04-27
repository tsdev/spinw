function mAtomList = matom(obj, varargin)
% generates all magnetic atoms in the unit cell
%
% Same as .atom, but only lists the magnetic atoms, which has non-zero
% spin. For speedup, call obj.matom(true), then it reads saved atomic
% positions from obj.
%

if (nargin>1) && ~isempty(obj.matomstore)
    
    mAtomList = obj.matomstore;
else
    
    atomList = obj.atom;
    
    mAtomList.r   = atomList.r(:,atomList.mag==1);
    mAtomList.idx = atomList.idx(:,atomList.mag==1);
    mAtomList.S   = obj.unit_cell.S(mAtomList.idx);
    
    obj.matomstore = mAtomList;
end

end % .matom