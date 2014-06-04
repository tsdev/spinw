function addg(obj, matrixIdx, varargin)
% assigns g-tensor to magnetic ions
%
% ADDG(obj, matrixIdx, {atomTypeIdx}, {atomIdx})
%
% Input:
%
% matrixIdx     Either an integer, that selects the matrix
%               obj.matrix.mat(:,:,matrixIdx), or a string identical to one
%               of the previously defined matrix labels, stored in
%               obj.matrix.label. Maximum value is nJ.
% atomTypeIdx   A vector that contains integers, the index of the magnetic
%               atoms in obj.unit_cell, with all symmetry equivalent atoms.
%               Maximum value is nAtom, if undefined g-tensor is assigned
%               to all magnetic atoms. Optional.
%  atomIdx      A vector that contains indices selecting some of the
%               symmetry equivalent atoms. Maximum value is the number of
%               symmetry equivalent atoms generated. If crystal symmetry is
%               not 0, atomIdx is not allowed, since the g-tensor for
%               equivalent atoms will be calculated using the symmetry
%               operators of the space group. Optional.
%
% Output:
%
% The function adds extra entries in the 'single_ion.g' field of the obj sw
% object.
%
% Example:
%
% ...
% cryst.addmatrix('label','g1','value',diag([1.8 1.8 2.1]))
% cryst.gencoupling
% cryst.addg('g1')
%
% This will add the 'g1' diagonal matrix to all magnetic atoms as
% anisotropic g-tensor.
%
% See also SW, SW.ADDCOUPLING, SW.ADDANISO, SW.ADDMATRIX.
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