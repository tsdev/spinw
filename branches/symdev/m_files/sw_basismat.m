function [M, asym] = sw_basismat(symOp, r, tol)
% [M, asym] = SW_BASISMAT(symOp, r, tol) determines the allowed matrix
% elements compatible with a given point group symmetry. The matrix can
% describe exchange interaction or single ion anisotropy.
%
% Input:
% symOp     Generators of the point group symmetry, dimensions are
%           [3 3 nSym]. Each symOp(:,:,ii) matrix defines a rotation.
% r         Distance vector between the two interacting atoms. For
%           anisotropy r=0, dimensions are [3 1].
% tol       Tolerance, optional. Default value is 1e-5.
%
% Output:
%
% M         Matrices, that span out the vector space of the symmetry
%           allowed matrices, dimensions are [3 3 nM]. Any matrix is
%           allowed that can be expressed as a linear combination of the
%           symmetry allowed matrices.
% asym      Logical vector, for each 3x3 matrix in M, tells whether is is
%           antisymmetry, dimensions are [1 nM].
%

if nargin == 0
    help sw_basismat;
    return;
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
        ordR = sw_rotorder(R);
        parR = 1;
    end
    if ordR == 10
        error('sw:sw_basismat:WrongOp','Not a valid point group generator!');
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
    
    M0.V = [dep(M0.V,M.V) dep(M.V,M0.V)];
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

rM = arrayfun(@(idx)rank([M0.V V0(:,idx)]),1:9);
% the nice part
rM = (rM==rank(M0.V));
Vnice = V0(:,rM);

% add rest of the vectors that cannot be expressed nicely :)
for ii = 1:size(M0.V,2)
    addV = indep(Vnice,M0.V(:,ii));
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

% reshape matrix
M = reshape(M,3,3,[]);

% selects the symmetric and antisymmetric matrices
asym = false(1,size(M,3));
for ii = 1:size(M,3)
    asym(ii) = sum(sum((M(:,:,ii)-M(:,:,ii)').^2)) > tol^2*9;
end

% sort assimetric matrices to the end
M = cat(3,M(:,:,~asym),M(:,:,asym));

% round M to the first 12th digit for better plotting 
% (slightly larger than the numerical error)
M = round(M*1e12)/1e12;

% cumprod(2*ones(1,9))
asym = [asym(~asym) asym(asym)];

end

function v = dep(M,v)
% returns column vectors of v that can be produced from the column vectors
% of M
%

rv = arrayfun(@(idx)rank([M v(:,idx)]),1:size(v,2));
rv = rv == rank(M);
v = v(:,rv);

end

function v = indep(M,v)
% returns column vectors of v that cannot be produced from the column
% vectors of M
%

rv = arrayfun(@(idx)rank([M v(:,idx)]),1:size(v,2));
rv = rv > rank(M);
v = v(:,rv);

end