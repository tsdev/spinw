function unit_cell_info = unit_cell(obj, idx)
% returns information on atoms.
%
% unit_cell_info = UNIT_CELL(obj, idx) returns information on symmetry
% inequivalent atoms. If idx is defined, it selects certain atoms, if it is
% undefined, information on all atoms are returned.
%
% Example:
% crystal1.unit_cell = unit_cell(crystal1,[1 3]);
% deletes the third atom of crystal1 sw object.
%

if nargin == 1
    idx = 1:numel(obj.unit_cell.S);
end

unit_cell_info.r     = obj.unit_cell.r(:,idx);
unit_cell_info.S     = obj.unit_cell.S(1,idx);
unit_cell_info.label = obj.unit_cell.label(1,idx);
unit_cell_info.color = obj.unit_cell.color(:,idx);

end