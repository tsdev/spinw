function M = sw_basismat(S, r)
% M = SW_BASISMAT(S, r) determines the allowed matrix elements compatible
% with a given point group symmetry. The matrix can describe exchange
% interaction or single ion anisotropy.
%
% Input:
% S         Generators of the point group symmetry, dimensions are
%           [3 3 nSym]. Each S(:,:,ii) matrix defines a rotation.
% r         Distance vector between the two interacting atoms. For
%           anisotropy r=0, dimensions are [1 3].
%
% Output:
%
% M         Matrices, that span out the vector space of the symmetry
%           allowed matrices, dimensions are [3 3 nM]. Any matrix is
%           allowed that can be expressed as a linear combination of the
%           symmetry allowed matrices.
%

nSym = size(S,3);

M0.V = zeros(9,0);
% normalize the r vector
if norm(r)>0
    r = r/norm(r);
end

for ii = 1:nSym
    % selected rotation operator
    R = S(:,:,ii);
    % order of the rotation operator
    ordR = 1;
    if norm(r)>0
        % find the number of rotation that overlay the two interacting atoms
        % (for anisotropy the 'interacting' atoms are identical)
        while (abs(abs(r'*(R^ordR)*r)-1)>1e-5) && (ordR<7)
            ordR = ordR + 1;
        end
    else
        while (norm(R^ordR-eye(3))>1e-5) && (ordR<7)
            ordR = ordR + 1;
        end
    end
    if ordR == 4
        error('sw:sw_basismat:WrongOp','Not a valid point group generator!');
    end
    % select the proper order that overlays the two atoms
    R = R^ordR;
    % parity of R in respect of the two interacting atoms
    parR = sign(r'*(R^ordR)*r);
    % solve the R*M-M*R=0 matrix equation valid for symmetric matrices
    [~, D, MS] = svd(kron(R,R)-eye(9));
    MS = MS(:,abs(diag(D))<1e-5);
    M.S = reshape(MS,3,3,[]);
    if parR==-1
        % solve the equation for antisymmetric matrices
        [~, D, MA] = svd(kron(R,R)+eye(9));
        MA = MA(:,abs(diag(D))<1e-5);
        M.A = reshape(MA,3,3,[]);
        M.A = M.A-permute(M.A,[2 1 3]);
    else
        M.A = M.S-permute(M.S,[2 1 3]);
    end
    M.S = M.S+permute(M.S,[2 1 3]);
    idx = [];
    M.V = cat(3,M.S,M.A);
    % remove zero vectors
    for jj = 1:size(M.V,3)
        if norm(M.V(:,:,jj))<1e-5
            idx = [idx jj]; %#ok<AGROW>
        end
    end
    M.V(:,:,idx) = [];
    M.V = reshape(M.V,9,[]);
    if ii>1
        % keep only those new vectors that are linearly dependent from the
        % eigenfunctions of the previous symmetry
        rM = arrayfun(@(idx)rank([M0.V M.V(:,idx)]),1:size(M.V,2));
        rM = rM==rank(M0.V);
        M0.V = M.V(:,rM);
    else
        M0.V = M.V;
    end
end

% make nice form 6 symmetric + 3 assymetric
%V0 = eye(9);
V0 = zeros(9);
% symmetric
V0(1,1) = 1; V0(5,2)=1; V0(9,3)=1; V0([2 4],4)=1; V0([3 7],5)=1; V0([6 8],6)=1;
% antisymmetric
V0([6 8],7)=[-1 1]; V0([3 7],8)=[1 -1]; V0([2 4],9)=[-1 1];

rM = arrayfun(@(idx)rank([M0.V V0(:,idx)]),1:9);
% the nice part
rM = (rM==rank(M0.V));
Vnice = V0(:,rM);
% add rest of the vectors that cannot be expressed nicely :)
rM = arrayfun(@(idx)rank([Vnice M0.V(:,idx)]),1:size(M0.V,2));
% the nice part
rM = (rank(Vnice)<rM);
M0.V = [Vnice M0.V(:,rM)];

% reshape matrix
M = reshape(M0.V,3,3,[]);
% normalize the maximum elements to one
nM  = arrayfun(@(idx)norm(M(:,:,idx)),1:size(M,3));
M = bsxfun(@rdivide,M,permute(nM,[1 3 2]));

end