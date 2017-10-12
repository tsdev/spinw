function [aList, SSext] = sw_extendlattice(nExt, aList, SS)
% creates superlattice
% 
% ### Syntax
% 
% `aList = sw_extendlattice(nExt,aList)`
% 
% `[aList, SSext] = sw_extendlattice(nExt,aList,SS)`
%
% ### Description
% 
% `aList = sw_extendlattice(nExt,aList)` creates a superlattice
% and calculates all atomic positions within the new superlattice by
% tiling it with the original cell.
%
% `[aList, SSext] = sw_extendlattice(nExt,aList,SS)` also calculates the
% bond matrix for the supercell by properly including all internal bonds
% and bonds between atoms in different supercells.
% 
% ### Input Arguments
% 
% `nExt`
% : Size of the supercell in units of the original cell in a row vector
%   with 3 elements.
% 
% `aList`
% : List of the atoms, produced by [spinw.matom].
% 
% `SS`
% : Interactions matrices in the unit cell. Struct where each field
%   contains an interaction matrix.
% 
% ### Output Arguments
% 
% `aList`
% : Parameters of the magnetic atoms in a struct with the following fields:
%   * `RRext` Positions of magnetic atoms in lattice units of the supercell stored in a matrix with dimensions of $[3\times n_{magExt}]$.
%   * `Sext`  Spin length of the magnetic atoms in a row vector with $n_{magExt}$ number of elements.
%
% `SSext`
% : Interaction matrix in the extended unit cell, struct type.
%   In the struct every field is a matrix. Every column of the
%   matrices describes a single bond, the following fields are generally
%   defined:
% 	* `iso`     Isotropic exchange interactions.
% 	* `ani`     Anisotropic exchange interations.
% 	* `dm`      Dzyaloshinsky-Moriya interaction terms.
% 	* `gen`     General $[3\times 3]$ matrix contains the exchange interaction.
% 
% ### See Also
% 
% [spinw.intmatrix]
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
    if ~isempty(SS.(fName))
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

end