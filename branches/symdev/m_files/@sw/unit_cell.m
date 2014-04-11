function unit_cell_info = unit_cell(obj, idx)
% returns information on atoms in the crystallographic unit cell.
%
% unit_cell_info = UNIT_CELL(obj, idx) returns information on symmetry
% inequivalent atoms. idx selects certain atoms, while UNIT_CELL(obj) or
% obj.UNIT_CELL returns information on all atoms. The selection can be also
% done according to the atom labels, in this case either a string of the
% label or cell of strings for several labels can be given.
%
% Sub fields are:
%   'r'         pasitions of the atoms in the unit cell, in a
%               3 x nAtom matrix, in lattice units
%   'S'         spin quantum number of the atoms, in a 1 x nAtom
%               vector, non-magnetic atoms have S=0
%   'label'     label of the atom, strings in a 1 x nAtom cell
%   'color'     color of the atom in 3 x nAtom matrix, where every
%               column is an 0-255 RGB color
%
% Example:
% crystal.unit_cell = unit_cell(crystal,[1 3]);
% keeps only the first and third symmetry inequivalent atoms in sw object.
%
% crystal.unit_cell = unit_cell(crystal,'O');
% keeps only the Oxygen atoms in crystal sw object.
%
% See also SW.ADDTWIN, SW.TWINQ, SW.UNIT_CELL.
%

if nargin == 1
    idx = 1:numel(obj.unit_cell.S);
else

% identify atoms by their labels
if ischar(idx)
    idx = {idx};
end

if iscell(idx)
    newIdx = [];
    for ii = 1:numel(idx)
        newIdx = [newIdx find(cellfun(@isempty,strfind(obj.unit_cell.label,idx{ii}))==false)]; %#ok<AGROW>
    end
    idx = newIdx;
end

end

unit_cell_info.r     = obj.unit_cell.r(:,idx);
unit_cell_info.S     = obj.unit_cell.S(1,idx);
unit_cell_info.label = obj.unit_cell.label(1,idx);
unit_cell_info.color = obj.unit_cell.color(:,idx);

end