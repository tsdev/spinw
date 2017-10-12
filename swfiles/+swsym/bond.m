function [genCp, ugenCp] = bond(r, bv, bond, symOp, tol)
% generates all symmetry equivalent bonds
%
% [genCp ugenCp] = SWSYM.BOND(r, bv, bond, symOp, {tol})
%
% Input:
%
% r         Positions of the magnetic atoms in the unit cell in a matrix
%           with dimensions of [3 nMagAtom], in lattice units.
% bv        Basis vectors of the lattice, used for error checking of the
%           calculation.
% bond      Vector, contains the starting bond with elements of 
%           [dl_x dl_y dl_z atom_1 atom_2], where dl is vector of lattice
%           translation between the two atom if they are not in the same
%           cell, atom_1 and atom_2 are indices of atoms in the list of
%           positions stored in r.
% symOp     Matrix, that contains the rotation and translation operators of
%           the space group with dimensions of [3 4 nOp].
% tol       Tolerance, optional. Default value is 1e-5.
%
% Output:
%
% genCp     Matrix, each column contains a coupling, the meaning of each
%           row are the same as the input coupling variable.
% ugenCp    Logical variable, true if the coupling is unique in the list of
%           generated couplings.
%
% See also SPINW.GENCOUPLING, SWSYM.OPERATOR, SWSYM.POSITION.
%

% TODO
tolDist = 1e-5;

if nargin == 0
    help swsym.bond
    return
end

if nargin < 4
    tol = 1e-5;
end

r1 = r(:,bond(4));
r2 = r(:,bond(5));
dl = bond(1:3);

% generate new atomic positions and translation vectors
r1new = permute(mmat(symOp(:,1:3,:),r1)+symOp(:,4,:),[1 3 2]);
r2new = permute(mmat(symOp(:,1:3,:),r2)+symOp(:,4,:),[1 3 2]);
dlnew = permute(mmat(symOp(:,1:3,:),dl),[1 3 2]) - cfloor(r1new,tol) + cfloor(r2new,tol);

% modulo to get atoms in the first unit cell
r1new = mod(r1new,1);
r2new = mod(r2new,1);

% throw away generated couplings with wrong distance
% determine the new indices in mAtom
[iNew, atom1] = isnewUC(r,r1new,tolDist);
if any(iNew)
    error('bond:SymProblem','The generated positions for atom1 are wrong!');
end
[iNew, atom2] = isnewUC(r,r2new,tolDist);
if any(iNew)
    error('bond:SymProblem','The generated positions for atom2 are wrong!');
end

dist = sqrt(sum((bv*(r(:,atom2)-r(:,atom1)+dlnew)).^2,1));
rightDist = abs(dist-dist(1)) < tol;
if ~all(rightDist)
    warning('bond:SymProblem','Symmetry generated couplings are dropped!');
end

genCp = [dlnew; atom1; atom2];
genCp = genCp(:,rightDist);

if nargout > 1
    % logical true if the coupling is unique in the list of generated couplings
    ugenCp = uniqueb(genCp);
end

end

function [isnew, symIdx] = isnewUC(A,B,tol)
% selects the new vectors from B within the first unit cell
%
% [isnew, symIdx] = isnewUC(A,B, tol)
%
% Dimensions of A and B have to be [3 nA] and [3 nB] respectively. A vector
% in B is considered new, if d(mod(vA-vB,1))<tol.
%
% Output:
%
% isnew     Vector of logical variables, true is the element of B differs
%           from all elements in A.
%

nA = size(A,2);
nB = size(B,2);

notequal = sum(sw_cmod(abs(repmat(permute(A,[2 3 1]),[1 nB 1]) - repmat(permute(B,[3 2 1]),[nA 1 1])),tol).^2,3) > tol;
isnew  = all(notequal,1);
idx    = 1:nB;
symIdx = max(bsxfun(@times,~notequal(:,idx(~isnew)),(1:size(notequal,1))'),[],1);

end

function r = sw_cmod(r, tol)
% modulo one with tolerance
% 
% ### Syntax
% 
% `r = sw_cmod(r, tol)`
% 
% ### Description
% 
% `r = sw_cmod(r, tol)` calculates modulo one with tolerance, numbers
% larger than $1-\epsilon$  $-\epsilon$.
% 
% ### See Also
% 
% [mod]
%

r = mod(r,1);

r(r > 1-tol) = r(r > 1-tol)-1;

end

function r = cfloor(r0, tol)
% floor for atomic positions
%
% floor(1-tol) == 1
%

r      = floor(r0);
idx    = abs(r0-r) > 1 - tol;
r(idx) = r(idx) + 1;

end

function uniqueB = uniqueb(bond)
% determines the unique bonds
%
% bond: [dl;atom1;atom2]
% two bonds are equivalent also when
% [dl;atom1;atom2] = [-dl;atom2;atom1]
%

nC  = size(bond,2);
c1  = permute(bond,[2 3 1]);
c2  = permute(bond,[3 2 1]);
nc1 = permute([-bond(1:3,:); bond([5 4],:)],[2 3 1]);

uniqueB = all(triu(any(bsxfun(@ne, c1,c2),3) & any(bsxfun(@ne,nc1,c2),3)) | tril(ones(nC)),1);

end