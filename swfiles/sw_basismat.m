function [M, asym] = sw_basismat(symOp, r, tol)
% determines allowed tensor components in a given point group symmetry
% 
% ### Syntax
% 
% `[m, asym] = sw_basismat(symop, r, tol)`
% 
% ### Description
% 
% `[m, asym] = sw_basismat(symop, r, tol)` determines the allowed tensor
% elements compatible with a given point group symmetry. The tensor can
% describe exchange interaction or single ion anisotropy. The function
% applies the symmetry invariance of the classical energy
% $\mathbf{S}_i\cdot \mathcal{M}\cdot \mathbf{S}_j$. Thus this symmetry
% analysis includes the transformation properties of spin operators as
% well.
% 
% ### Input Arguments
% 
% `symOp`
% : Generators of the point group symmetry, in a matrix with dimensions of
%   $[3\times 3\times n_{sym}]$ where each `symOp(:,:,ii)` matrix defines a rotation.
% 
% `r`
% : Distance column vector between the two interacting atoms. For
%   anisotropy $r=0$.
% 
% `tol`
% : Tolerance, optional, default value is $10^{-5}$.
% 
% ### Output Arguments
% 
% `M`
% : Matrices, that span out the vector space of the symmetry
%           allowed matrices, dimensions are $[3\times 3\times n_M]$. Any matrix is
%           allowed that can be expressed as a linear combination of the
%           symmetry allowed matrices.
%
% `asym`
% : Logical vector, for each $[3\times 3]$ matrix in $M$, tells whether it is
%           antisymmetric stored in a row vector with $n_M$ elements.
% 
% ### See Also
% 
% [spinw.getmatrix] \| [spinw.setmatrix]
%

if nargin == 0
    swhelp sw_basismat
    return
elseif nargin == 2
    tol = 1e-5;
end

nSym = size(symOp,3);

% make nice form 6 symmetric + 3 asymetric matrices as starting bases
V0 = zeros(9);
% symmetric
V0(1,1) = 1; V0(5,2)=1; V0(9,3)=1; V0([2 4],4)=1; V0([3 7],5)=1; V0([6 8],6)=1;
% antisymmetric
V0([6 8],7)=[-1 1]; V0([3 7],8)=[1 -1]; V0([2 4],9)=[-1 1];

%
M0.V = V0;
% normalize the r vector
if norm(r)>0
    r = r/norm(r);
    aniso = false;
else
    aniso = true;
end

for ii = 1:nSym
    % selected rotation operator
    R = symOp(:,:,ii);
    if ~aniso
        % find the number of rotation that overlay the two interacting atoms
        % (for anisotropy the 'interacting' atoms are identical)
        ordR = 1;
        while (abs(abs(r'*(R^ordR)*r)-1) > tol) && (ordR<10)
            ordR = ordR + 1;
        end
        % select the proper order that overlays the two atoms
        R = R^ordR;
        % parity of R in respect of the two interacting atoms
        parR = sign(r'*R*r);
        
    else
        % check that the rotation operator is valid
        ordR = swsym.oporder([R zeros(3,1)]);
        parR = 1;
    end
    if ordR == 10
        error('sw_basismat:WrongOp','Not a valid point group generator!');
    end
    % solve the R*M-M*R=0 matrix equation valid for symmetric matrices
    [~, D, MS] = svd(kron(R,R)-eye(9));
    M.S = reshape(MS(:,abs(diag(D)) < tol),3,3,[]);
    
    if parR == -1
        % solve the equation for antisymmetric matrices
        [~, D, MA] = svd(kron(R,R)+eye(9));
        M.A = reshape(MA(:,abs(diag(D)) < tol),3,3,[]);
        M.A = M.A-permute(M.A,[2 1 3]);
    else
        M.A = M.S-permute(M.S,[2 1 3]);
    end
    M.S = M.S+permute(M.S,[2 1 3]);
    M.V = reshape(cat(3,M.S,M.A),9,[]);
    
    %M0.V = [dep(M0.V,M.V) dep(M.V,M0.V)];
    M0.V = intersec(M0.V,M.V);
    
    % remove small matrices
    normM = arrayfun(@(idx)norm(M0.V(:,idx)),1:size(M0.V,2));
    M0.V(:,normM < tol) = [];
    
end

% separate symmetric and antisymmetric components
M0.V = reshape(M0.V,3,3,[]);
M0.V = reshape(cat(3,M0.V-permute(M0.V,[2 1 3]),M0.V+permute(M0.V,[2 1 3])),9,[]);

% remove small matrices
normM = arrayfun(@(idx)norm(M0.V(:,idx)),1:size(M0.V,2));
M0.V(:,normM < tol) = [];

rM = arrayfun(@(idx)rank([M0.V V0(:,idx)],tol),1:9);
% the nice part
rM = (rM==rank(M0.V,tol));
Vnice = V0(:,rM);

% add rest of the vectors that cannot be expressed nicely :)
for ii = 1:size(M0.V,2)
    addV = indep(Vnice,M0.V(:,ii),tol);
    if ~isempty(addV)
        addV = addV - sum(bsxfun(@times,sum(bsxfun(@times,addV,Vnice),1),Vnice),2);
        Vnice = [Vnice addV]; %#ok<AGROW>
    end
end

%M0.V = orth([Vnice M0.V]);
%M = M0.V;
M = Vnice;

% normalize the largest absolute value element to one
divM  = arrayfun(@(idx)M(find(abs(M(:,idx))==max(abs(M(:,idx))),1,'first'),idx),1:size(M,2));
M = bsxfun(@rdivide,M,divM);

% subtract vectors from each other to get nicer result
for ii = 1:(size(M,2)-1)
    for jj = (ii+1):size(M,2)
        factM = M(:,jj).*M(:,ii);
        if any(factM)
            factM = sum(factM)/sum(factM~=0);
        else
            factM = 0;
        end
         M(:,jj) = M(:,jj)-M(:,ii)*factM;
    end
end

% remove small matrices
normM = arrayfun(@(idx)norm(M(:,idx)),1:size(M,2));
M(:,normM < tol) = [];

% normalize the largest absolute value element to one
divM  = arrayfun(@(idx)M(find(abs(M(:,idx))==max(abs(M(:,idx))),1,'first'),idx),1:size(M,2));
M = bsxfun(@rdivide,M,divM);

% subtract vectors from each other to get nicer result backwards as well
M = M(:,end:-1:1);
for ii = 1:(size(M,2)-1)
    for jj = (ii+1):size(M,2)
        factM = M(:,jj).*M(:,ii);
        if any(factM)
            factM = sum(factM)/sum(factM~=0);
        else
            factM = 0;
        end
         M(:,jj) = M(:,jj)-M(:,ii)*factM;
    end
end

% normalize the largest absolute value element to one
divM  = arrayfun(@(idx)M(find(abs(M(:,idx))==max(abs(M(:,idx))),1,'first'),idx),1:size(M,2));
M = bsxfun(@rdivide,M,divM);

% reshape matrix
M = reshape(M,3,3,[]);

% selects the symmetric and antisymmetric matrices
asym = false(1,size(M,3));
for ii = 1:size(M,3)
    asym(ii) = sum(sum((M(:,:,ii)-M(:,:,ii)').^2)) > tol^2*9;
end

aIdx = find(asym);

% sort assimetrix matrices to Dx, Dy, Dz order
if any(asym)
    [~, aSort] = sort(squeeze(abs(M(2,3,asym)) + abs(M(1,3,asym))*10 + abs(M(1,2,asym))*100));
    M(:,:,aIdx) = M(:,:,aIdx(aSort));
end
% choose proper sign of the DM interactions
for ii = 1:numel(aIdx)
    signM = sign(M(2,3,aIdx(ii))*100 + M(3,1,aIdx(ii))*10 + M(1,2,aIdx(ii)));
    M(:,:,aIdx(ii)) = M(:,:,aIdx(ii)) * signM;
end



% sort assimetric matrices to the end
M = cat(3,M(:,:,~asym),M(:,:,asym));

% round M to the first 12th digit for better plotting
% (slightly larger than the numerical error)
M = round(M*1e12)/1e12;

% cumprod(2*ones(1,9))
asym = [asym(~asym) asym(asym)];

end

% function v = dep(M,v,tol)
% % returns column vectors of v that can be produced from the column vectors
% % of M
% %
% 
% rv = arrayfun(@(idx)rank([M v(:,idx)],tol),1:size(v,2));
% rv = rv == rank(M,tol);
% v = v(:,rv);
% 
% end

function I = intersec(U1,U2)
% returns the intersection of two vectorspaces U1 and U2. Each column of U1
% and U2 are a vector.

% N = basis for nullspace of [U1 U2]
N = null([U1 U2]); 
% I = basis for intersection of U1 and U2
I = U1*N(1:size(U1,2),:); 
% I = orthonormal and minimal size basis
I = orth(I); 

end


function v = indep(M,v, tol)
% returns column vectors of v that cannot be produced from the column
% vectors of M
%

rv = arrayfun(@(idx)rank([M v(:,idx)],tol),1:size(v,2));
rv = rv > rank(M,tol);
v = v(:,rv);

end