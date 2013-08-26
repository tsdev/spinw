function addg(obj, matrixIdx, varargin)
% assigns g-tensor to magnetic ions
%
% ADDG(obj, matrixIdx, {atomTypeIdx}, {atomIdx})
%
% matrixIdx         Either an integer, that selects one of the matrices in
%                   obj.matrix.mat, or a string that identical to one of
%                   the previously defined labels, stored in
%                   obj.matrix.label. Maximum value is nCoupling.
% atomTypeIdx       A vector that contains integers, that point to the atom
%                   in obj.unit_cell, then it selects atoms and
%                   symmetry equivalent ones with coordinates:
%                       unit_cell.r(:,atomIdx)
%                   Maximum value is nAtom, it is optional, if undefined
%                   g-tensor is assigned to all magnetic atoms. Optional.
%  atomIdx          A vector that contains indices of the symmetry
%                   equivalent atoms, selecting only a few of them. Maximum
%                   value is the number of symmetry equivalent atoms
%                   generated. If crystal symmetry is not P1, atomIdx is
%                   not allowed, since the g-tensor for equivalent
%                   atoms will be calculated using the symmetry operators
%                   of the space group. Optional.
%
% See also SW.ADDCOUPLING, SW.ADDANISO.
%

mAtom    = obj.matom;
nMagAtom = size(mAtom.idx,2);

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
        error('sw:addg:SymmetryProblem','atomIdx is not allowed when crystal symmetry is not P1!');
    end

end

addField = obj.single_ion;

if nargin > 2
    atomTypeIdx = varargin{1};
    
    if length(obj.single_ion.g) ~= nMagAtom
        addField.g = zeros(nMagAtom,1);
    end
    
    for ii = 1:length(atomTypeIdx)
        aTemp = addField.g(mAtom.idx == atomTypeIdx(ii));
        if nargin > 3
            aTemp(atomIdx) = matrixIdx;
        else
            aTemp = matrixIdx;
        end
        addField.g(mAtom.idx == atomTypeIdx(ii)) = aTemp;
    end
else
    addField.g = zeros(1,nMagAtom) + matrixIdx;
end

addField.g = int32(addField.g);
obj.single_ion = addField;

end