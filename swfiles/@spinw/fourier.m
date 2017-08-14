function res = fourier(obj,hkl,varargin)
% calculates the Fourier transformation of the Hamiltonian
%
% res = FOURIER(obj,hkl,'option1', value1 ...)
%
% The function calculates the following sum:
%       J(k) = sum_ij J_ij * exp(i*k*d_ij)
% The code is optimised for calculating the sum for large number of
% k-vectors and alternatively for a large number of d_ij. The single ion
% anisotropy is not included in the sum.
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
% extend        If true, the Fourier transform will be calculated on the
%               magnetic supercell, if false the crystallographic cell will
%               be considered. Default is true.
% isomode       Defines how Heisenberg/non-Heisenberg Hamiltonians are
%               treated. Can have the following values:
%                   'off'   Always output the (3x3) form of the
%                           Hamiltonian, (default).
%                   'auto'  If the Hamiltonian is Heisenberg, only output
%                           one of the diagonal values from the (3x3)
%                           matrices to reduce memory consumption.
% fid           Defines whether to provide text output. Default is defined
%               by the swpref.getpref('fid') command. The possible values:
%                   0       No text output is generated.
%                   1       Text output in the MATLAB Command Window.
%                   fid     File ID provided by the fopen() command, the
%                           output is written into the opened file stream.
%
% Output:
%
% res           Structure with the following fields:
%   ft          contains the Fourier transfor in a matrix with dimensions
%               [3,3,nMagExt,nMagExt,nHKL] or [1,1,nMagExt,nMagExt,nHKL]
%               for Heisenberg and non-Heisenberg Hamiltonians respectively
%               (if isomode is 'auto'). Here nMagExt is the number of
%               magnetic atoms in the magnetic cell and nHKL is the number
%               of reciprocal space points.
%   hkl         Matrix with the given reciprocal space points stored in a
%               matrix with dimensions [3,nHKL].
%   isiso       True is the output is in Heisenberg mode, when the ft
%               matrix has dimensions of [1,1,nMagExt,nMagExt,nHKL],
%               otherwise is false.
%
% See also SPINW.OPTMAGK.
%

% TODO: test for magnetic supercell

inpForm.fname  = {'extend' 'fid' 'isomode' 'sublat'};
inpForm.defval = {true     nan   'off'     []      };
inpForm.size   = {[1 1]    [1 1] [1 -1]    [1 -2]  };
inpForm.soft   = {false    false false     true    };

param = sw_readparam(inpForm, varargin{:});

switch param.isomode
    case 'auto'
        isomode = true;
    case 'off'
        isomode = false;
    otherwise
        error('spinw:WrongInput','The given isomode option is invalid!')
end

if isnan(param.fid)
    % Print output into the following file
    fid = swpref.getpref('fid',true);
else
    fid = param.fid;
end

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

% calculate symbolic Fourier transformation if obj is in symbolic mode
if obj.symbolic
    if numel(hkl) == 3
        hkl = sym(hkl);
    end
    
    if ~isa(hkl,'sym')
        inpForm.fname  = {'fitmode'};
        inpForm.defval = {false    };
        inpForm.size   = {[1 1]    };
        param0 = sw_readparam(inpForm, varargin{:});
        
        if ~param0.fitmode
            fprintf0(obj.fileid,['No symbolic hkl value was given, Fourier'...
                ' transformation for general Q (h,k,l) will be calculated!\n']);
        end
        res = obj.fouriersym(varargin{:});
    else
        res = obj.fouriersym(varargin{:},'hkl',hkl);
    end
    return
end

% generate exchange couplings
%[SS, ~, RR] = obj.intmatrix('fitmode',true,'extend',param.extend,'conjugate',true);
[SS, ~, RR] = obj.intmatrix(struct('fitmode',true,'extend',param.extend,'conjugate',true),'noCheck');

% is there Heisenberg exchange only
matidx = unique(obj.coupling.mat_idx(:)');
matidx(matidx==0) = [];
isIso = all(sw_mattype(obj.matrix.mat(:,:,matidx))==1);

fprintf0(fid,'Calculating the Fourier transformation of the spin Hamiltonian...\n')
if isIso
    fprintf0(fid,'Heisenberg model...\n')
else
    fprintf0(fid,'Non-Heisenberg model...\n')
end
% list of magnetic atoms in the unit cell
matom = obj.matom;
% number of Q points
nHkl    = size(hkl,2);
if param.extend
    nExt = double(obj.mag_str.nExt);
else
    nExt = 1;
end
% number of magnetic atoms in the crystallographic unit cell
nMag    = numel(matom.idx);
% number of bonds
nBond   = size(SS.all,2);

% interacting atom1
atom1 = mod(SS.all(4,:)-1,nMag)+1;
% interacting atom2
atom2 = mod(SS.all(5,:)-1,nMag)+1;

% magnetic couplings, 3x3xnJ
%JJ = cat(3,reshape(SS.all(6:14,:),3,3,[]),SI.aniso);

if isIso
    % keep only a single value for Heisenberg interactions
    JJ = SS.all(6,:);
else
    JJ = SS.all(6:14,:);
end

% include the spin values in the exchange energy: J_ij*S_i*S_j
JJ = bsxfun(@times,JJ,matom.S(atom1).*matom.S(atom2));

% distance vector between interacting atoms in l.u.
dR = SS.all(1:3,:)+RR(:,atom2)-RR(:,atom1);

% exponents, dimension: 1 x nHkl x nBond
ExpF = permute(exp(1i*2*pi*sum(bsxfun(@times,hkl',permute(dR,[3 1 2])),2)),[2 1 3]);

% J*exp(ikr), dimensions 9 x nBond x nHkl (1 x nBond x nHkl for Heisenberg)
Jexp = permute(bsxfun(@times,permute(JJ,[1 3 2]),ExpF),[1 3 2]);

% reduce the number of sublattices
if ~isempty(param.sublat)
    atom1 = param.sublat(atom1);
    atom2 = param.sublat(atom2);
    nMag  = max(param.sublat);
    fprintf0(fid,'Remapping magnetic atoms into a new set of sublattices...\n');
end

% number of magnetic atoms in the magnetic supercell
nMagExt = prod(nExt)*nMag;

% optimize for large number of hkl, is not faster for more than ~10
% bonds% if nHkl>nBond
%     ft = zeros(3,3,nMagExt,nMagExt,nHkl);
%     if isIso
%         for ii = 1:nBond
%             ft(1,1,atom1(ii),atom2(ii),:) = ft(1,1,atom1(ii),atom2(ii),:)+permute(Jexp(1,ii,:),[1 2 4 5 3]);
%         end
%         ft(2,2,:,:,:) = ft(1,1,:,:,:);
%         ft(3,3,:,:,:) = ft(1,1,:,:,:);
%     else
%         for ii = 1:nBond
%             ft(:,:,atom1(ii),atom2(ii),:) = reshape(Jexp(:,ii,:),3,3,[]);
%         end
%     end
% end

% optimize for large number of bonds
if ~isIso
    % general, non-Heisenberg model
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
    % Dimensions: 9 x nMagExt x nMagExt x nHkl (1 x nMagExt x nMagExt x nHkl for Heisenberg)
    %                 atom1     atom2
    % the result is always real, because we use the conjugate of the bonds as
    % well
    ft = real(accumarray(idxAll,Jexp(:)',[9 nMagExt nMagExt nHkl]));
    % reshape into 3 x 3 x nMagExt x nMagExt x nHkl
    ft = reshape(ft,[3 3 nMagExt nMagExt nHkl]);
else
    % Heisenberg model
    % J^\alpha^\beta component indices is always 1
    % ...
    % atom1 indices
    idx2 = repmat(atom1,[1 1 nHkl]);
    % atom2 indices
    idx3 = repmat(atom2,[1 1 nHkl]);
    % hkl indices
    idx4 = repmat(permute(1:nHkl,[1 3 2]),[1 nBond 1]);
    % indices all
    idxAll = [idx2(:) idx3(:) idx4(:)];
    
    % Exchange energy: A_ij = S_i*S_j*J_ij*exp(ik*2pi*(R_i-R_j))
    % Summed up on all bonds
    % Dimensions: 9 x nMagExt x nMagExt x nHkl (1 x nMagExt x nMagExt x nHkl for Heisenberg)
    %                 atom1     atom2
    % the result is always real, because we use the conjugate of the bonds as
    % well
    ft = real(accumarray(idxAll,Jexp(:)',[nMagExt nMagExt nHkl]));
    % reshape into 3 x 3 x nMagExt x nMagExt x nHkl
    ft = reshape(ft,[1 1 nMagExt nMagExt nHkl]);
    if ~isomode
        ft(3,3,:,:,:) = ft(1,1,:,:,:);
        ft(2,2,:,:,:) = ft(1,1,:,:,:);
    end
end

% save results in a struct
res.ft    = ft;
res.hkl   = hkl;
% Heisenberg output
res.isiso = size(ft,1)==1;

fprintf0(fid,'Calculation finished.\n');

end