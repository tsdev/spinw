function res = fourier(obj,hkl)
% calculates the Fourier transformation of the Hamiltonian
%
% res = FOURIER(obj,hkl)
%
%

% TODO: test for magnetic supercell

% for linear scans create the Q line(s)
if nargin > 1
    if iscell(hkl)
        hkl = sw_qscan(hkl);
    elseif numel(hkl)==3
        hkl = hkl(:);
    end
else
    hkl = [];
end

% generate exchange couplings
[SS, ~, RR] = obj.intmatrix('fitmode',2,'extend',true,'conjugate',false);

% list of magnetic atoms in the unit cell
matom = obj.matom;
% number of Q points
nHkl    = size(hkl,2);
% number of magnetic atoms in the magnetic supercell
nMagExt = prod(obj.mag_str.N_ext)*numel(matom.idx);
% number of bonds
nBond   = size(SS.all,2);

% interacting atom1
atom1 = SS.all(4,:);
% interacting atom1
atom2 = SS.all(5,:);

% exchange energy: J_ij*S_i*S_j
JJ = bsxfun(@times,SS.all(6:14,:),matom.S(atom1).*matom.S(atom2));

% distance vector between interacting atoms in l.u.
dR = SS.all(1:3,:)+RR(:,SS.all(5,:))-RR(:,SS.all(4,:));

% exponents, dimension: 1 x nHkl x nBond
ExpF = permute(sum(exp(1i*2*pi*bsxfun(@times,hkl',permute(dR,[3 1 2]))),2),[2 1 3]);

% J*exp(ikr), dimensions 9 x nHkl x nBond
Jexp = permute(bsxfun(@times,permute(JJ,[1 3 2]),ExpF),[1 3 2]);

% J^\alpha^\beta component indices in numbers 1:9
idx1 = repmat((1:9)',[1 nBond nHkl]);
% atom1 indices
idx2 = repmat(atom1,[9 1 nHkl]);
% atom2 indices
idx3 = repmat(atom2,[9 1 nHkl]);
% hkl indices
idx4 = repmat(permute(1:nHkl,[1 3 2]),[9 nBond 1]);
% indices all
idxAll = [idx1(:) idx2(:) idx3(:) idx4(:)];

% Exchange energy: A_ij = S_i*S_j*J_ij*exp(ik*2pi*(R_i-R_j))
% Summed up on all bonds
% Dimensions: 9 x nMagExt x nMagExt x nHkl
%                 atom1     atom2
ft = accumarray(idxAll,Jexp(:)',[9 nMagExt nMagExt nHkl]);

% reshape into 3 x 3 x nMagExt x nMagExt x nHkl
ft = reshape(ft,[3 3 nMagExt nMagExt nHkl]);

% add complex conjugate
ft = ft + conj(permute(ft,[1 2 4 3 5]));
%ft = real(ft);
%ft = ft + conj(permute(ft,[1 3 2 4]));

% save results in a struct
res.ft  = ft;
res.hkl = hkl;

end