function res = fouriersym(obj,varargin)
% calculates the Fourier transformation of a symbolic Hamiltonian
%
% res = FOURIER(obj, 'option1', value1 ...)
%
% Input:
%
% obj           Input structure, spinw class object.
%
% Options:
%
% hkl           Symbolic definition of q vector. Default is the general Q
%               point:
%                   hkl = [sym('h') sym('k') sym('l')]
%
%
%
% See also SPINW.FOURIER.
%

% TODO: test for magnetic supercell

hkl0 = [sym('h','real'); sym('k','real'); sym('l','real')];

title0 = 'Symbolic Fourier transformation of the Hamiltonian';

inpForm.fname  = {'tol' 'hkl'  'title'};
inpForm.defval = {1e-4   hkl0  title0 };
inpForm.size   = {[1 1] [3 1]  [1 -1] };

param = sw_readparam(inpForm, varargin{:});

% symbolic hkl 
hkl = param.hkl;

% generate exchange couplings
[SS, ~, RR] = obj.intmatrix('fitmode',2,'extend',true,'conjugate',true);

% list of magnetic atoms in the unit cell
matom = obj.matom;
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
dR = SS.all(1:3,:) + (RR(:,SS.all(5,:))-RR(:,SS.all(4,:)));

% exponents, dimension: 1 x nBond
ExpF = permute(exp(1i*2*pi*sum(bsxfun(@times,hkl',permute(dR,[3 1 2])),2)),[1 3 2]);

% J*exp(ikr), dimensions 3 x 3 x nBond
Jexp = reshape(bsxfun(@times,JJ,ExpF),3,3,nBond);

% Exchange energy: A_ij = S_i*S_j*J_ij*exp(ik*2pi*(R_i-R_j))
% Summed up on all bonds
% Dimensions: 3 x 3 x nMagExt x nMagExt
%                     atom1     atom2

ft = sym(zeros(3,3,nMagExt,nMagExt));
for ii = 1:nBond
    ft(:,:,atom1(ii),atom2(ii)) = ft(:,:,atom1(ii),atom2(ii)) + Jexp(:,:,ii);
end
    
% save results in a struct
res.ft  = ft;
res.hkl = hkl;

end