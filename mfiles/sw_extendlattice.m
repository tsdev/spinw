function [aList, SSext] = sw_extendlattice(nExt, aList, SS)
% creates superlattice
%
% [aList, SSext] = SW_EXTENDLATTICE(nExt, aList, {SS})
%
% It creates a superlattice and all redefines all given bond for the larger
% superlattice.
%
% Input:
%
% nExt          Number of unit cell extensions, dimensions are [1 3].
% aList         List of the atoms, produced by sw.matom.
% SS            Interactions matrices in the unit cell, optional.
%
% Output:
%
% aList         Parameters of the magnetic atoms.
% aList.RRext   Positions of magnetic atoms, assuming an extended unit
%               cell, dimensions are [3 nMagExt].
% aList.Sext    Spin length of the magnetic atoms, dimensions are
%               [1 nMagExt].
%
% SSext         Interaction matrix in the extended unit cell, struct type.
%               In the struct every field is a matrix. Every column of the
%               matrices describes a single interaction.
% SSext.iso     Isotropic exchange interactions.
% SSext.ani     Anisotropic exchange interations.
% SSext.dm      Dzyaloshinsky-Moriya interaction terms.
% SSext.gen     General 3x3 matrix contains the exchange interaction.
%
% See also SW.INTMATRIX.
%

if nargin == 0
    help sw_extendlattice
    return
end

nAtom    = size(aList.r,2);
nCell    = prod(nExt);
nExt     = nExt(:);
nExt1    = nExt-1;

% generate cell indices in the supercell
[cIdx{1:3}] = ndgrid(0:nExt1(1),0:nExt1(2),0:nExt1(3));
cIdx        = cat(4,cIdx{:});

% generate spin quantum numbers for the extendd unit cell
aList.Sext = repmat(aList.S,[1 nCell]);
% generate atomic positions in the extended unit cell
aList.RRext = bsxfun(@plus,bsxfun(@rdivide,cIdx,permute(nExt,[2:4 1])),permute(bsxfun(@rdivide,aList.r,nExt),[3 4 5 1 2]));
aList.RRext = reshape(permute(aList.RRext,[4 5 1:3]),3,[]);
% generate the indices for each atom in the extended unit cell
if isfield(aList,'idx')
    aList.idxext = repmat(aList.idx,[1 nCell]);
end

if nargin < 3
    SSext = struct;
    return
end

% loop over all field names of SS
fNameV = fields(SS);
SSext = struct;

% additional atom index for bonds [1 1 nCell*nAtom]
addIdx = permute((0:(nCell-1))*nAtom,[1 3 2]);

for ii = 1:numel(fNameV)
    fName = fNameV{ii};
    SSext.(fName) = repmat(SS.(fName),[1 nCell]);
    % first atom index within the uspercell
    SSext.(fName)(4,:) = reshape(bsxfun(@plus,addIdx,SS.(fName)(4,:)),1,[]);
    % end of bond vector still in original cell dimensions
    bVect = reshape(permute(bsxfun(@plus,cIdx,permute(SS.(fName)(1:3,:),[3:5 1 2])),[4 5 1:3]),3,[]);
    % normalize bond vector to supercell dimensions
    SSext.(fName)(1:3,:) = floor(bsxfun(@rdivide,bVect,nExt));
    % indices are between (0:nCell-1)*nAtom
    SSext.(fName)(5,:) = sum(bsxfun(@times,bVect-bsxfun(@times,SSext.(fName)(1:3,:),nExt),[1;nExt(1);prod(nExt(1:2))]),1)*nAtom+SSext.(fName)(5,:);
end

end