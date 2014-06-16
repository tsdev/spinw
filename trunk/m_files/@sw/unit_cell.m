function unit_cell_info = unit_cell(obj, idx)
% returns information on atoms in the crystallographic unit cell
%
% unit_cell_info = UNIT_CELL(obj, idx)
%
% The function returns information on symmetry inequivalent atoms. 
%
% Input:
%
% obj       sw class object.
% idx       Selects certain atoms. If undefined UNIT_CELL(obj) or
%           obj.UNIT_CELL returns information on all atoms. The selection
%           can be also done according to the atom labels, in this case
%           either a string of the label or cell of strings for several
%           labels can be given.
%
% Output:
%
% 'unit_cell_info' is a tructure with the following fields:
% r         Positions of the atoms in the unit cell, in a matrix with
%           dimensions of [3 nAtom], in lattice units.
% S         Spin quantum number of the atoms, in a [1 nAtom] horizontal
%           vector, non-magnetic atoms have S=0.
% label     Label of the atom, strings in a cell with dimensions of 
%           [1 nAtom].
% color     Color of the atom in a matrix with dimensions of [3 nAtom],
%           where every column is an 0-255 RGB color code.
%
% Example:
%
% ...
% cryst.unit_cell = unit_cell(cryst,[1 3]);
%
% The example keeps only the first and third symmetry inequivalent atoms in
% cryst object.
%
% ...
% cryst.unit_cell = unit_cell(cryst,'O');
%
% The example keeps only the Oxygen atoms in cryst object.
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