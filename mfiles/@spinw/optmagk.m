function [hkl, Eopt, n, stat] = optmagk(obj,varargin)
% determines the magnetic propagation vector
%
% res = OPTMAGK(obj,hkl,'option1', value1 ...)
%
% Input:
%
% obj           Input structure, spinw class object.
% hkl           Defines the Q points where the Fourier transform is
%               calculated, in reciprocal lattice units, size is [3 nHkl].
%               Q can be also defined by several linear scan in reciprocal
%               space. In this case hkl is cell type, where each element of
%               the cell defines a point in Q space. Linear scans are
%               assumed between consecutive points. Also the number of Q
%               points can be specified as a last element, it is 100 by
%               defaults. For example: hkl = {[0 0 0] [1 0 0]  50}, defines
%               a scan along (h,0,0) from 0 to 1 and 50 Q points are
%               calculated along the scan.
%
%               For symbolic calculation at a general reciprocal space
%               point use sym class input. For example to calculate the
%               spectrum along (h,0,0): hkl = [sym('h') 0 0]. To do
%               calculation at a specific point do for example sym([0 1
%               0]), to calculate the spectrum at (0,1,0).
%
% Options:
%
% fitmode       Speedup (for fitting mode only), default is false.
%
%
% See also SPINW.TISZA.
%

% TODO: test for magnetic supercell

%inpForm.fname  = {'fitmode' };
%inpForm.defval = {false     };
%inpForm.size   = {[1 1]     };
%
%param = sw_readparam(inpForm, varargin{:});



% calculate symbolic Fourier transformation if obj is in symbolic mode
if obj.symbolic
    warning('spinw:optmagk:NoSymbolic','The function does not work in symbolic mode!');
    return
end


% generate exchange couplings
[SS, SI, RR] = obj.intmatrix('fitmode',2,'extend',false,'conjugate',true);

% list of magnetic atoms in the unit cell
matom = obj.matom;
% number of magnetic atoms in the unit cell
nMagAtom = numel(matom.idx);
% number of bonds
nBond   = size(SS.all,2);

% interacting atom1
%atom1 = SS.all(4,:);
atom1 = mod(SS.all(4,:)-1,numel(matom.S))+1;
% interacting atom1
%atom2 = SS.all(5,:);
atom2 = mod(SS.all(5,:)-1,numel(matom.S))+1;

% exchange energy: J_ij*S_i*S_j
JJ = bsxfun(@times,SS.all(6:14,:),matom.S(atom1).*matom.S(atom2));

% distance vector between interacting atoms in l.u.
dR = SS.all(1:3,:)+RR(:,atom2)-RR(:,atom1);

% number of Q points
nHkl    = 1;

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

% optimise the energy
[hkl, ~, stat] = ndbase.pso([],@optfun,[1/4 1/4 1/4],'lb',[0 0 0],'ub',[1/2 1/2 1/2],varargin{:});

[Eopt, V] = optfun(hkl);

n = cross(real(V(1:3)),imag(V(1:3)));
n = n/norm(n);
n = n(:)';

    function [E, V] = optfun(p)
        %
        
        % exponents, dimension: 1 x nHkl x nBond
        ExpF = permute(exp(1i*2*pi*sum(bsxfun(@times,p(:)',permute(dR,[3 1 2])),2)),[2 1 3]);
        
        % J*exp(ikr), dimensions 9 x nBond x nHkl
        Jexp = permute(bsxfun(@times,permute(JJ,[1 3 2]),ExpF),[1 3 2]);
        
        % Exchange energy: A_ij = S_i*S_j*J_ij*exp(ik*2pi*(R_i-R_j))
        % Summed up on all bonds
        % Dimensions: 9 x nMagExt x nMagExt x nHkl
        %                 atom1     atom2
        ft = accumarray(idxAll,Jexp(:)',[9 nMagAtom nMagAtom nHkl]);
        
        % reshape into 3 x 3 x nMagExt x nMagExt x nHkl
        E = min(eig(reshape(permute(reshape(ft,[3 3 nMagAtom nMagAtom nHkl]),[1 3 2 4 5]),3*nMagAtom,3*nMagAtom,[])));
        
        if nargout > 1
            [V,~] = eig(reshape(permute(reshape(ft,[3 3 nMagAtom nMagAtom nHkl]),[1 3 2 4 5]),3*nMagAtom,3*nMagAtom,[]));
            V = V(:,1);
        end
    end

end