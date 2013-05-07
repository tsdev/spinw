function spectra = spinwave(obj, hkl, varargin)
% calculates dynamical spin-spin correlation function using linear spin wave theory
%
% spectra = SPINWAVE(obj, k, 'option1', value1 ...)
%
% Spin wave dispersion and spin-spin correlation function is calculated at
% the reciprocal space points k. The function can deal with arbitrary
% commensurate magnetic structure and magnetic interactions as well as
% single ion anisotropy and magnetic field.
%
% Input:
%
% obj           Input structure, sw class object.
% hkl           Defines the Q points where the spectra is calculated, in
%               reciprocal lattice units, size is [3 nHkl]. Q can be also
%               defined by several linear scan in reciprocal space. In this
%               case hkl is cell type, where each element of the cell
%               defines a point in Q space. Linear scans are assumed
%               between consecutive points. Also the number of Q points can
%               be specified as a last element, it is 100 by defaults. For
%               example: hkl = {[0 0 0] [1 0 0]  50}, defines a scan along
%               (h,0,0) from 0 to 1 and 50 Q points are calculated along
%               the scan.
%
% Options:
%
% formfact      TODO !!!!!!!!!!!!!!!!!!!!!!!!
%               Whether to include the magnetic form factor in the swConv
%               convoluted spectra. If true the form factor based on the
%               name of the magnetic atoms are used to read form factor
%               data from formfactor.dat file.
% fitmode       Speedup (for fitting mode only), default is false. In
%               fitmode the eigenvalue sorting is skipped, see
%               <a href="matlab:doc eigenshuffle">eigenshuffle</a>.
% notwin        If true, the spectra of the twins won't be calculated.
%               Default is false.
%
% Output:
%
% spectra is a structure, with the following fields: 
% omega         Calculated spin wave dispersion, dimensins are 
%               [2*nMagExt nHkl], where nMagExt is the number of magnetic
%               atoms in the extended unit cell.
% Sab           Dynamical structure factor, dimensins are 
%               [3 3 2*nMagExt nHkl]. Each (:,:,i,j) submatrix contains the
%               9 correlation functions: Sxx, Sxy, Sxz, etc.
%
% If several domains exist in the sample, omega and Sab are packaged into a
% cell, that contains nTwin number of matrices.
%
% hkl           Contains the input Q values, dimensins are [3 nHkl].
% hklA          Same Q values, but in reciproc Angstrom units in the
%               lab coordinate system, dimensins are [3 nHkl].
% ff            Magnetic form factor for all magnetic ions calculated
%               for hklA momentum transfer values, dimensins are 
%               [nMagExt nHkl].
% obj           The copy of the input obj.
%
% See also SW, SW_NEUTRON, SW.SWINC, SW_POL, SW.POWSPEC, SW.OPTMAGSTR.
%

% help when executed without argument
if nargin==1
    help sw.spinwave
    return
end


% for linear scans create the Q line(s)
if iscell(hkl)
    hkl = sw_qscan(hkl);
end

inpForm.fname  = {'fitmode' 'notwin'};
inpForm.defval = {false     false   };
inpForm.size   = {[1 1]     [1 1]   };

param = sw_readparam(inpForm, varargin{:});

nExt    = double(obj.mag_str.N_ext);
spectra = struct;
nHkl    = size(hkl,2);

% Calculates momentum transfer in A^-1 units.
hklA = 2*pi*(hkl'/obj.basisvector)';

% define Q scans for the twins
nTwin = size(obj.twin.vol,2);
if param.notwin
    nTwin = 1;
end

nHkl0 = nHkl;
if nTwin>1
    % In the abc coordinate system of the selected twin the scan is
    % rotated opposite direction to rotC.
    hkl  = cell2mat(obj.twinq(hkl));
    nHkl = nHkl*nTwin;
end

% Create the interaction matrix and atomic positions in the extended
% magnetic unit cell.
if param.fitmode
    [SS, SI, RR] = obj.intmatrix('fitmode',2);
else
    [SS, SI, RR] = obj.intmatrix;
end

% Introduce the opposite couplings.
% (i-->j) and (j-->i)
% transpose the JJ matrix as well [1 2 3 4 5 6 7 8 9] --> [6 9 12 7 10 13 8 11 14]
SS.new         = [SS.all(1:3,:)   -SS.all(1:3,:)  ];
SS.new(4:5,:)  = [SS.all([4 5],:)  SS.all([5 4],:)];
SS.new(6:14,:) = [SS.all(6:14,:)   SS.all([6 9 12 7 10 13 8 11 14],:) ]/2;
SS.all         = SS.new;

% Converts wavevctor list into the extended unit cell
hklExt = bsxfun(@times,hkl,nExt')*2*pi;

% Calculates parameters eta and zed.
if isempty(obj.mag_str.S)
    error('sw:sw_spinwave:NoMagneticStr','No magnetic structure defined in obj!');
end

M0  = obj.mag_str.S;
S0  = sqrt(sum(M0.^2,1));
nMagExt = size(M0,2);

% Local (e1,e2,e3) coordinate system fixed to the moments.
% e3 || Si
e3 = bsxfun(@rdivide,M0,S0);
% e2 = Si x [1,0,0], if Si || [1,0,0] --> e2 = [0,0,1]
e2  = [zeros(1,nMagExt); e3(3,:); -e3(2,:)];
%e2(3,~any(e2)) = 1;
e2(3,~any(abs(e2)>1e-10)) = 1;
e2  = bsxfun(@rdivide,e2,sqrt(sum(e2.^2,1)));
% e1 = e2 x e3
e1  = cross(e2,e3);

% Defines eta and zed.
zed = e1 + 1i*e2;
eta = e3;

dR    = [SS.all(1:3,:) zeros(3,nMagExt)];
atom1 = [SS.all(4,:)   1:nMagExt];
atom2 = [SS.all(5,:)   1:nMagExt];
JJ    = cat(3,reshape(SS.all(6:end,:),3,3,[]),SI.aniso);
nCoupling = size(JJ,3);

zedL = repmat(permute(zed(:,atom1),[1 3 2]),[1 3 1]);
zedR = repmat(permute(zed(:,atom2),[3 1 2]),[3 1 1]);

etaL = repmat(permute(eta(:,atom1),[1 3 2]),[1 3 1]);
etaR = repmat(permute(eta(:,atom2),[3 1 2]),[3 1 1]);

SiSj = sqrt(S0(atom1).*S0(atom2));

% Creates temporary values for calculating matrix elements.
AD  =  shiftdim(sum(sum(etaL.*JJ.*etaR,2),1),1);
A20 = -S0(atom2).*AD;
D20 = -S0(atom1).*AD;
BC0 =  SiSj.*shiftdim(sum(sum(zedL.*JJ.*     zedR ,2),1),1);
AD0 =  SiSj.*shiftdim(sum(sum(zedL.*JJ.*conj(zedR),2),1),1);
MF  =  repmat(obj.unit.gamma*SI.field*eta,[1 2]);

% Creates the serial indices for every matrix element in ham matrix.
idxA1 = [atom1'         atom2'         ];
idxA2 = [atom1'         atom1'         ];
idxB  = [atom1'         atom2'+nMagExt ];
idxC  = [atom1'+nMagExt atom2'         ];
idxD1 = idxA1+nMagExt;
idxD2 = [atom2'+nMagExt atom2'+nMagExt ];
idxMF = [(1:2*nMagExt)' (1:2*nMagExt)' ];

% Creates the matrix of exponential factors nCoupling x nHkl size.
% Extends dR into 3 x 3 x nCoupling x nHkl
ExpF = exp(1i*permute(sum(repmat(dR,[1 1 nHkl]).*repmat(permute(hklExt,[1 3 2]),[1 nCoupling 1]),1),[2 3 1]))';

% Creates the matrix elements containing zed.
A1 = bsxfun(@times,     AD0 ,ExpF);
B  = bsxfun(@times,     BC0 ,ExpF);
C  = bsxfun(@times,conj(BC0),ExpF);
D1 = bsxfun(@times,conj(AD0),ExpF);

% Store all indices
idxAll = [idxA1; idxB; idxC; idxD1];
% Store all matrix elements
ABCD   = [A1     B     C     D1   ];

% Stores the matrix elements in ham.
idx3   = repmat(1:nHkl,[4*nCoupling 1]);
idxAll = [repmat(idxAll,[nHkl 1]) idx3(:)];
idxAll = idxAll(:,[2 1 3]);

ABCD   = ABCD';

ham = accumarray(idxAll,ABCD(:),[2*nMagExt 2*nMagExt nHkl]);

ham = ham + repmat(accumarray([idxA2; idxD2],2*[A20 D20],[1 1]*2*nMagExt),[1 1 nHkl]);

if any(SI.field)
    ham = ham + repmat(accumarray(idxMF,MF,[1 1]*2*nMagExt),[1 1 nHkl]);
end

ham = (ham + conj(permute(ham,[2 1 3])))/2;

g = diag([ones(nMagExt,1); -ones(nMagExt,1)]);
% Dispersion of all spin wave modes, size: 2*nMagExt x nHkl.
omega = zeros(2*nMagExt,nHkl);

% All the matrix calculations are according to White's paper

gham = 0*ham;
for ii = 1:nHkl
    gham(:,:,ii) = g*ham(:,:,ii);
end

if param.fitmode
    % Without eigenshuffle
    V = zeros(2*nMagExt,2*nMagExt,nHkl);
    D = zeros(2*nMagExt,nHkl);
    
    for ii = 1:nHkl
        [V(:,:,ii), Dtemp] = eig(gham(:,:,ii));
        D(:,ii)     = diag(Dtemp);
    end
else
    % With eigenshuffle
    [V,D] = eigenshuffle(gham);
end

for ii = 1:nHkl
    omega(:,ii) = diag(g).*D(:,ii);
    M           = diag(g*V(:,:,ii)'*g*V(:,:,ii));
    V(:,:,ii)   = V(:,:,ii)*diag(sqrt(1./M));
end

% Calculates correlation functions.
VExtR = repmat(permute(V      ,[4 5 1 2 3]),[3 3 1 1 1]);
VExtL = repmat(permute(conj(V),[4 5 2 1 3]),[3 3 1 1 1]);

% Introduces the exp(-ikR) exponential factor.
ExpF =  exp(-1i*sum(repmat(permute(hklExt,[1 3 2]),[1 nMagExt 1]).*repmat(RR,[1 1 nHkl]),1));
% Includes the sqrt(Si/2) prefactor.
ExpF = ExpF.*repmat(sqrt(S0/2),[1 1 nHkl]);

ExpFL = repmat(permute(ExpF,[1 4 5 2 3]),[3 3 2*nMagExt 2]);
ExpFR = conj(repmat(permute(ExpF,[1 4 2 5 3]),[3 3 2 2*nMagExt]));

zeda = repmat(permute([zed conj(zed)],[1 3 4 2]),[1 3 2*nMagExt 1         nHkl]);
zedb = repmat(permute([conj(zed) zed],[3 1 2])  ,[3 1 1         2*nMagExt nHkl]);

% Dynamical for factor from S^alpha^beta(k) correlation function.
% Sab(alpha,beta,iMode,iHkl), size: 3 x 3 x 2*nMagExt x nHkl.
% Normalizes the intensity to single unit cell.
Sab = squeeze(sum(zeda.*ExpFL.*VExtL,4)).*squeeze(sum(zedb.*ExpFR.*VExtR,3))/prod(nExt);

if nTwin>1
    % Rotate the calculated correlation function into the twin coordinate
    % system using rotC
    SabAll = cell(1,nTwin);
    for ii = 1:nTwin
        % select the ii-th twin from the Q points
        idx    = (1:nHkl0) + (ii-1)*nHkl0;
        % select correlation function of twin ii
        SabT   = Sab(:,:,:,idx);
        % size of the correlation function matrix
        sSabT  = size(SabT);
        % convert the matrix into cell of 3x3 matrices
        SabT   = reshape(SabT,3,3,[]);
        % select the rotation matrix of twin ii
        rotC   = obj.twin.rotc(:,:,ii);
        % rotate correlation function using arrayfun
        SabRot = arrayfun(@(idx)(rotC*SabT(:,:,idx)*(rotC')),1:size(SabT,3),'UniformOutput',false);
        SabRot = cat(3,SabRot{:});
        % resize back the correlation matrix
        SabAll{ii} = reshape(SabRot,sSabT);
    end
    Sab = SabAll;
    omega = mat2cell(omega,2*nMagExt,repmat(nHkl0,[1 nTwin]));
end

% Creates output structure with the calculated values.
spectra.omega  = omega;
spectra.Sab    = Sab;
spectra.hkl    = hkl(:,1:nHkl0);
spectra.hklA   = hklA;

if ~param.fitmode
    spectra.ff     = ones(nMagExt,nHkl0);
    spectra.obj    = copy(obj);
end

end