function unit_cell_info = unitcell(obj, idx)
% returns unit cell data
% 
% ### Syntax
% 
% `cellInfo = unitcell(obj, idx)`
% 
% ### Description
% 
% `cellInfo = unitcell(obj, idx)` returns information on symmetry
% inequivalent atoms and allowing to subselect certain atoms using the
% `idx` index vector.
% 
% ### Examples
% 
% The example keeps only the first and third symmetry inequivalent atoms in
% `cryst` object.
% ```
% cryst.unit_cell = unitcell(cryst,[1 3]);
% ```
% The example keeps only the atoms with labels `'O'` (Oxygen) atoms in
% `cryst` object.
% ```
% cryst.unit_cell = unitcell(cryst,'O');
% ```
% 
% ### Input Arguments
% 
% `obj`
% : [spinw] object.
% 
% `idx`
% : Selects certain atoms. If undefined `unit_cell(obj)` or
%      `obj.unit_cell` returns information on all atoms. The selection
%      can be also done according to the atom labels, in this case
%      either a string of the label or cell of strings for several
%      labels can be given.
% 
% ### Output Arguments
% 
% `cellInfo`
% : Structure that contains all the fields of [spinw.unit_cell].
% 
% ### See Also
% 
% [spinw.addtwin] \| [spinw.twinq] \| [spinw.unit_cell]
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

% copy all fields
%fName = {'r' 'S' 'label' 'color' 'ox' 'occ' 'b' 'ff' 'A' 'Z' 'biso'};
fName = fieldnames(obj.unit_cell);
for ii = 1:numel(fName)
    unit_cell_info.(fName{ii}) = obj.unit_cell.(fName{ii})(:,idx);
end
unit_cell_info.ff = obj.unit_cell.ff(:,:,idx);

end