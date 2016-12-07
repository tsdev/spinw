function [genCp, ugenCp] = sw_gensymcoupling(obj, coupling, sym, tol)
% generates all equivalent couplings, using space group symmetry
%
% [genCp ugenCp] = SW_GENSYMCOUPLING(obj, coupling, sym, {tol})
%
% Input:
%
% obj       spinw class object.
% coupling  Vector, contains the original coupling. Dimensions are [5 1].
%           Elements are [dlx dly dlz atom1 atom2 idx dist].
% sym       Cell, that contains the rotation and translation operators of
%           the space group: {rotOp, trOp}.
% tol       Tolerance, optional. Default value is 1e-5.
%
% Output:
%
% genCp     Matrix, each column contains a coupling, the meaning of each
%           row are the same as the input coupling variable.
% ugenCp    Logical variable, true if the coupling is unique in the list of
%           generated couplings.
%

% TODO
tolDist = 1e-5;

if nargin == 0
    help sw_gensymcoupling
    return
end

if nargin < 4
    tol = 1e-5;
end

mAtom = obj.matom;

symOp = sym{1};
symTr = sym{2};

r1 = mAtom.r(:,coupling(4));
r2 = mAtom.r(:,coupling(5));
dl = coupling(1:3);

% generate new atomic positions and translation vectors
r1new = permute(mmat(symOp,r1),[1 3 2]) + symTr;
r2new = permute(mmat(symOp,r2),[1 3 2]) + symTr;
%dlnew = permute(mmat(symOp,dl),[1 3 2]) - ceil(r1new) + floor(r2new) + 1;
dlnew = permute(mmat(symOp,dl),[1 3 2]) - cfloor(r1new,tol) + cfloor(r2new,tol);
% modulo to get atoms in the first unit cell
r1new = mod(r1new,1);
r2new = mod(r2new,1);
% throw away generated couplings with wrong distance
% determine the new indices in mAtom
[iNew, atom1] = isnewUC(mAtom.r,r1new,tolDist);
if any(iNew)
    error('Sym error: generated positions for atom1 are wrong!');
end
[iNew, atom2] = isnewUC(mAtom.r,r2new,tolDist);
if any(iNew)
    error('Sym error: generated positions for atom2 are wrong!');
end
dist = sqrt(sum((obj.basisvector*(mAtom.r(:,atom2)-mAtom.r(:,atom1)+dlnew)).^2,1));

rightDist = abs(dist-dist(1)) < tol;
if ~all(rightDist)
    warning('Symmetry generated couplings are dropped!');
end

genCp = [dlnew; atom1; atom2];
genCp = genCp(:,rightDist);

if nargout > 1
    % logical true if the coupling is unique in the list of generated couplings
    ugenCp = uniquec(genCp);
end

end

function [isnew, symIdx] = isnewUC(A,B, tol)
% [isnew, symIdx] = isnewUC(A,B, tol)
% selects the new vectors from B within the first unit cell. Dimensions of
% A and B have to be [3 nA] and [3 nB] respectively.
% A vector in B is considered new, if d(mod(vA-vB,1))<tol.
%
% Output:
%
% isnew     Vector of logical variables, true is the element of B differs
%           from all elements in A.
%

nA = size(A,2);
nB = size(B,2);

notequal = sum(sw_cmod(abs(repmat(permute(A,[2 3 1]),[1 nB 1]) - repmat(permute(B,[3 2 1]),[nA 1 1])),tol).^2,3) > tol;

isnew = all(notequal,1);

idx = 1:nB;

%symIdx = arrayfun(@(idx)find(~notequal(:,idx),1,'first'),idx(~isnew));
% faster solution
symIdx = max(bsxfun(@times,~notequal(:,idx(~isnew)),(1:size(notequal,1))'),[],1);

end

function r = cfloor(r0, tol)
% floor for atomic positions:
% floor(1-tol) == 1
%

r = floor(r0);

idx = abs(r0-r) > 1 - tol;

r(idx) = r(idx) + 1;

end

function uniqueC = uniquec(coupling)
% determines the unique couplings
% coupling: [dl;atom1;atom2]
% two couplings are equivalent also when
% [dl;atom1;atom2] = [-dl;atom2;atom1]
%

nC = size(coupling,2);
c1 = permute(coupling,[2 3 1]);
c2 = permute(coupling,[3 2 1]);
nc1 = permute([-coupling(1:3,:); coupling([5 4],:)],[2 3 1]);

uniqueC = all(triu(any(bsxfun(@ne, c1,c2),3) & any(bsxfun(@ne,nc1,c2),3)) | tril(ones(nC)),1);

end