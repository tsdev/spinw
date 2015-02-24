function addaniso(obj, matrixIdx, varargin)
% assigns anisotropy matrices to magnetic ions
%
% ADDANISO(obj, matrixIdx, {atomTypeIdx}, {atomIdx})
%
% Input:
%
% matrixIdx     Either an integer, that selects the matrix
%               obj.matrix.mat(:,:,matrixIdx), or a string identical to one
%               of the previously defined matrix labels, stored in
%               obj.matrix.label. Maximum value is nJ.
% atomTypeIdx   String or cell of strings that select magnetic atoms by
%               their label. Also can be a vector that contains integers,
%               the index of the magnetic atoms in obj.unit_cell, with all
%               symmetry equivalent atoms. Maximum value is nAtom, if
%               undefined anisotropy is assigned to all magnetic atoms.
%               Optional.
%  atomIdx      A vector that contains indices selecting some of the
%               symmetry equivalent atoms. Maximum value is the number of
%               symmetry equivalent atoms generated. If crystal symmetry is
%               not 0, atomIdx is not allowed, since the anisotropy matrix
%               for equivalent atoms will be calculated using the symmetry
%               operators of the space group. Optional.
%
% Output:
%
% The function adds extra entries in the 'single_ion.aniso' field of the
% obj sw object.
%
% Example:
%
% ...
% cryst.addmatrix('label','A1','value',diag([-0.1 -0.1 0]))
% cryst.gencoupling
% cryst.addaniso('A1')
%
% This will add the 'A1' diagonal matrix to all magnetic atoms as
% anisotropy (easy XY plane anisotropy).
%
% See also SW, SW.ADDCOUPLING, SW.ADDG, SW.ADDMATRIX.
%

mAtom    = obj.matom;
nMagAtom = size(mAtom.idx,2);

if nMagAtom==0
    error('sw:addaniso:NoMagAtom','There is no magnetic atom in the unit cell with S>0!');
end

cTemp = -1;
if ischar(matrixIdx)
    for ii = length(obj.matrix.label):-1:1
        if strcmp(obj.matrix.label{ii},matrixIdx)
            cTemp = ii;
        end
    end
    matrixIdx = cTemp;
    
end

if matrixIdx == -1
    error('sw:addaniso:WrongCouplingTypeIdx','Input matrix does not exists!');
end

if nargin > 3
    atomIdx = varargin{2};
    if obj.lattice.sym > 1
        error('sw:addaniso:SymmetryProblem','atomIdx is not allowed when crystal symmetry is not P1!');
    end

end

addField = obj.single_ion;

if nargin > 2
    atomTypeIdx = varargin{1};
    
    if length(obj.single_ion.aniso) ~= nMagAtom
        addField.aniso = zeros(nMagAtom,1);
    end
    
    % select atoms by label
    if ischar(atomTypeIdx)
        atomTypeIdx = {atomTypeIdx};
    end
    
    if iscell(atomTypeIdx)
        % loop over all atom labels
        isSelectedAtom = zeros(1,obj.natom);
        for ii = 1:numel(atomTypeIdx)
            % atom index
            sAIdx = find(strcmp(obj.unit_cell.label,atomTypeIdx{ii}));
            if numel(sAIdx) == 0
                error('sw:addaniso:WrongString','Atom label doesn''t exists!')
            elseif numel(sAIdx) > 1
                error('sw:addaniso:WrongString','Multiple atoms have the same label!')
            end
            isSelectedAtom(sAIdx) = true;
        end
        atomTypeIdx = find(isSelectedAtom);
    end

    for ii = 1:length(atomTypeIdx)
        aTemp = addField.aniso(mAtom.idx == atomTypeIdx(ii));
        if nargin > 3
            aTemp(atomIdx) = matrixIdx;
        else
            aTemp = matrixIdx;
        end
        addField.aniso(mAtom.idx == atomTypeIdx(ii)) = aTemp;
    end
else
    addField.aniso = zeros(1,nMagAtom)+matrixIdx;
end

addField.aniso = int32(addField.aniso);
obj.single_ion = addField;

end