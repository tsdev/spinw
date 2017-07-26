function unit_cell_info = unitcell(obj, idx)
% returns information on atoms in the crystallographic unit cell
%
% unit_cell_info = UNITCELL(obj, idx)
%
% The function returns information on symmetry inequivalent atoms. 
%
% Input:
%
% obj       spinw class object.
% idx       Selects certain atoms. If undefined UNIT_CELL(obj) or
%           obj.UNIT_CELL returns information on all atoms. The selection
%           can be also done according to the atom labels, in this case
%           either a string of the label or cell of strings for several
%           labels can be given.
%
% Output:
%
% 'unit_cell_info' is a tructure with that contains all the fields of
% unit_cell.
%
% Example:
%
% ...
% cryst.unit_cell = unitcell(cryst,[1 3]);
%
% The example keeps only the first and third symmetry inequivalent atoms in
% cryst object.
%
% ...
% cryst.unit_cell = unitcell(cryst,'O');
%
% The example keeps only the Oxygen atoms in cryst object.
%
% See also SPINW.ADDTWIN, SPINW.TWINQ, SPINW.UNIT_CELL.
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