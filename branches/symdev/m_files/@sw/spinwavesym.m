function spectra = spinwavesym(obj, varargin)
% calculates dynamical spin-spin correlation function using linear spin wave theory
%
% spectra = SPINWAVESYM(obj, k, 'option1', value1 ...)
%
% Spin wave dispersion and spin-spin correlation function is calculated at
% the reciprocal space points k. The function can deal with arbitrary
% magnetic structure and magnetic interactions as well as single ion
% anisotropy and magnetic field.
%
% If the magnetic ordering wavevector is non-integer, the dispersion is
% calculated using a coordinate system rotating from cell to cell. In this
% case the Hamiltonian has to fulfill this extra rotational symmetry.
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
% fitmode       Speedup (for fitting mode only), default is false.
% notwin        If true, the spectra of the twins won't be calculated.
%               Default is false.
% fid           For text output of the calculation:
%               0   No text output.
%               1   Output written onto the Command Window. (default)
%               fid Output written into a text file opened with the
%                   fid = fopen(path) command.
% tol           Tolerance of the incommensurability of the magnetic
%               ordering wavevector. Deviations from integer values of the
%               ordering wavevector smaller than the tolerance are
%               considered to be commensurate. Default value is 1e-4.
% omega_tol     Tolerance on the energy difference of degenerate modes when
%               diagonalising the quadratic form, default is 1e-5.
% hermit        Method for matrix diagonalization:
%                   true        Method of J.H.P. Colpa, Physica 93A (1978)
%                               327.
%                   false       Method of R.M. White, PR 139 (1965) A450.
%               Colpa: the grand dynamical matrix is converted into another
%                      Hermitian matrix, that can will give the real
%                      eigenvalues.
%               White: the non-Hermitian g*H matrix will be diagonalised,
%                      that is strictly speaking not the right solution
%                      method.
%               Advise:
%               Always use Colpa's method, except when small imaginary
%               eigenvalues are expected. In this case only White's method
%               work. The solution in this case is wrong, however by
%               examining the eigenvalues it can give a hint where the
%               problem is.
%               Default is true.
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
% incomm        Whether the spectra calculated is incommensurate or not.
% obj           The copy of the input obj.
%
% See also SW, SW.SPINWAVE, SW_NEUTRON, SW_POL, SW.POWSPEC, SW.OPTMAGSTR.
%

% help when executed without argument
if nargin==1
    help sw.spinwave
    return
end

inpForm.fname  = {'fitmode' 'notwin' 'sortMode' 'optmem' 'fid' 'tol' 'omega_tol' 'hermit'};
inpForm.defval = {false     false    true      true     1     1e-4   1e-5        true    };
inpForm.size   = {[1 1]     [1 1]    [1 1]      [1 1]    [1 1] [1 1]  [1 1]      [1 1]   };

param = sw_readparam(inpForm, varargin{:});

if param.fitmode
    param.sortMode = false;
end

% seize of the extended magnetic unit cell
nExt    = double(obj.mag_str.N_ext);
% magnetic ordering wavevector in the extended magnetic unit cell
km = obj.mag_str.k.*nExt;
% whether the structure is incommensurate
incomm = any(abs(km-round(km)) > param.tol);

% symbolic wavevectors
syms h k l real
hkl = [h; k; l];

fid = param.fid;

spectra = struct;

% Create the interaction matrix and atomic positions in the extended
% magnetic unit cell.
[SS, SI, RR] = obj.intmatrix('plotmode',true,'extend',true);


% Introduce the opposite couplings.
% (i-->j) and (j-->i)
% transpose the JJ matrix as well [1 2 3 4 5 6 7 8 9] --> [6 9 12 7 10 13 8 11 14]
SS.new         = [SS.all(1:3,:)   -SS.all(1:3,:)  ];
SS.new(4:5,:)  = [SS.all([4 5],:)  SS.all([5 4],:)];
SS.new(6:14,:) = [SS.all(6:14,:)   SS.all([6 9 12 7 10 13 8 11 14],:) ]/2;

% ???
idxM = [SS.all(15,:) SS.all(15,:)];

SS.all         = SS.new;

% Converts wavevctor list into the extended unit cell
hklExt  = hkl.*nExt'*2*pi;

% Calculates parameters eta and zed.
if isempty(obj.mag_str.S)
    error('sw:spinwave:NoMagneticStr','No magnetic structure defined in obj!');
end

M0 = obj.mag_str.S;
S0 = sqrt(sum(M0.^2,1));
% normal to rotation of the magnetic moments
n  = obj.mag_str.n;
nMagExt = size(M0,2);

if fid ~= 0
    if incomm
        fprintf(fid,'Calculating SYMBOLIC INCOMMENSURATE spin wave spectra (nMagExt = %d)...\n',nMagExt);
    else
        fprintf(fid,'Calculating SYMBOLIC COMMENSURATE spin wave spectra (nMagExt = %d)...\n',nMagExt);
    end
end

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
idxM  = [idxM repmat(obj.single_ion.aniso,1,prod(nMagExt))];
% magnetic couplings, 3x3xnJ
JJ = cat(3,reshape(SS.all(6:14,:),3,3,[]),SI.aniso);

% remove zero anisotropy matrices
anyIdx = squeeze(sum(sum(abs(JJ),1),2))>0;
dR    = dR(:,anyIdx);
atom1 = atom1(1,anyIdx);
atom2 = atom2(1,anyIdx);
JJ    = JJ(:,:,anyIdx);
idxM  = idxM(1,anyIdx);

% create symbolic variables from matrix labels
nMat = size(obj.matrix.mat,3);

for ii = 1:nMat
    JS(ii) = sym(obj.matrix.label{ii},'real'); %#ok<AGROW>
end

% normalize JJ to the maximum value
JJ = JJ./repmat(max(max(abs(JJ))),[3 3 1]);


% create symbolic matrices
JJ = JJ.*repmat(permute(JS(idxM),[3 1 2]),[3 3 1]);

if incomm
    % transform JJ due to the incommensurate wavevector
    [~, K] = sw_rot(n,km*dR*2*pi);
    % multiply JJ with K matrices for every interaction
    % and symmetrising JJ for the rotating basis
    JJ = (mmat(JJ,K)+mmat(K,JJ))/2;
end

zedL = repmat(permute(zed(:,atom1),[1 3 2]),[1 3 1]);
zedR = repmat(permute(zed(:,atom2),[3 1 2]),[3 1 1]);

etaL = repmat(permute(eta(:,atom1),[1 3 2]),[1 3 1]);
etaR = repmat(permute(eta(:,atom2),[3 1 2]),[3 1 1]);

SiSj = sqrt(S0(atom1).*S0(atom2));

% Creates temporary values for calculating matrix elements.
A0 =  SiSj.*shiftdim(symsum2(symsum2(zedL.*JJ.*conj(zedR),2),1),1);
B0 =  SiSj.*shiftdim(symsum2(symsum2(zedL.*JJ.*     zedR ,2),1),1);
C0 =  shiftdim(symsum2(symsum2(etaL.*JJ.*etaR,2),1),1);
C1 = -2*S0(atom2).*C0;
C2 = -2*S0(atom1).*C0;

% Magnetic field
MF =  repmat(obj.unit.muB*SI.field*eta,[1 2]);

% Creates the serial indices for every matrix element in ham matrix.
idxA1 = [atom1'         atom2'         ];
idxA2 = idxA1+nMagExt;

idxC1 = [atom1'         atom1'         ];
idxC2 = [atom2'+nMagExt atom2'+nMagExt ];

idxB1 = [atom1'         atom2'+nMagExt ];
% transpose of idxB
idxB2 = idxB1(:,[2 1]);

idxMF = [(1:2*nMagExt)' (1:2*nMagExt)' ];

% Creates the matrix of exponential factors nCoupling x 1 size.
ExpF = exp(1i*hklExt'*dR);

% Creates the matrix elements containing zed.
A1 = A0.*ExpF;
A2 = conj(A0).*ExpF;
B1 = B0.*ExpF;
B2 = conj(B1);

% Store all indices
idxAll = [idxA1; idxB1; idxB2; idxA2; idxC1; idxC2; idxMF];
% Store all matrix elements
ABCD   = [A1     B1     B2     A2     C1     C2     MF   ];

ham = sym(zeros(2*nMagExt));

for ii = 1:size(idxAll,1)
    ham(idxAll(ii,1),idxAll(ii,2)) = ham(idxAll(ii,1),idxAll(ii,2)) + ABCD(ii);
end

ham = simplify((ham + conj(permute(ham,[2 1 3])))/2);

g = sym(diag([ones(nMagExt,1); -ones(nMagExt,1)]));

if param.hermit
    % All the matrix calculations are according to Colpa's paper
    % J.H.P. Colpa, Physica 93A (1978) 327-353
    
%         posDef = feval(symengine,'linalg::isPosDef',ham);
%         
%         if ~posDef
%             ham = ham + sym(eye(2*nMagExt)*param.omega_tol);
%             posDef = feval(symengine,'linalg::isPosDef',ham);
%             warning('sw:spinwavesym',['To make ham positive definite, a '...
%                 'small omega_tol added to the diagonal!'])
%         end
        
        K = simplify(feval(symengine,'linalg::factorCholesky',ham,'NoCheck'));
        
%         if ~posDef
%             warning('sw:spinwavesym:PositiveDefiniteHamiltonian',...
%                 ['Hamiltonian matrix is not positive definite, probably'...
%                 ' the magnetic structure is wrong!']);
%         end
        
        K2 = K*g*K';
        K2 = 1/2*(K2+K2');
        
        % Hermitian K2 will give orthogonal eigenvectors
        [U, omega] = eig(K2);
        
        % the inverse of the para-unitary transformation V
        V = inv(K)*U*diag(sqrt(g*omega)); %#ok<MINV>
        
        omega = diag(omega);

else
    % All the matrix calculations are according to White's paper
    % R.M. White, et al., Physical Review 139, A450?A454 (1965)
    
    [V, D] = eig(g*ham);
    
    % multiplication with g removed to get negative and positive
    % energies as well
    omega = diag(D);
    M = diag(g*V'*g*V);
    V = V*diag(sqrt(1./M));

end


% Calculates correlation functions.
% V right
VExtR = repmat(permute(V      ,[4 5 1 2 3]),[3 3 1 1 1]);
% V left: conjgate transpose of V
VExtL = conj(permute(VExtR,[1 2 4 3 5]));
%VExtL = repmat(permute(conj(V),[4 5 2 1 3]),[3 3 1 1 1]);

% Introduces the exp(-ikR) exponential factor.
ExpF =  exp(-1i*sum(repmat(permute(hklExt0MEM,[1 3 2]),[1 nMagExt 1]).*repmat(RR,[1 1 nHklMEM]),1));
% Includes the sqrt(Si/2) prefactor.
ExpF = ExpF.*repmat(sqrt(S0/2),[1 1 nHklMEM]);

ExpFL =      repmat(permute(ExpF,[1 4 5 2 3]),[3 3 2*nMagExt         2]);
% conj transpose of ExpFL
ExpFR = conj(permute(ExpFL,[1 2 4 3 5]));
%ExpFR = conj(repmat(permute(ExpF,[1 4 2 5 3]),[3 3 2         2*nMagExt]));

zeda = repmat(permute([zed conj(zed)],[1 3 4 2]),[1 3 2*nMagExt 1         nHklMEM]);
% conj transpose of zeda
zedb = conj(permute(zeda,[2 1 4 3 5]));
%zedb = repmat(permute([conj(zed) zed],[3 1 2 4]),[3 1 1         2*nMagExt nHklMEM]);

% Dynamical for factor from S^alpha^beta(k) correlation function.
% Sab(alpha,beta,iMode,iHkl), size: 3 x 3 x 2*nMagExt x nHkl.
% Normalizes the intensity to single unit cell.
Sab = cat(4,Sab,squeeze(sum(zeda.*ExpFL.*VExtL,4)).*squeeze(sum(zedb.*ExpFR.*VExtR,3))/prod(nExt));

if fid ~= 0
    fprintf(fid,'Calculation finished.\n');
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% END MEMORY MANAGEMENT LOOP
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if incomm
    % resize matrices due to the incommensurability (k-km,k,k+km) multiplicity
    kmIdx = repmat(sort(repmat([1 2 3],1,nHkl0/3)),1,nTwin);
    % Rodrigues' rotation formula.
    nx  = [0 -n(3) n(2); n(3) 0 -n(1); -n(2) n(1) 0];
    nxn = n'*n;
    K1 = 1/2*(eye(3) - nxn - 1i*nx);
    K2 = nxn;
    
    % symmetrise Sab by averaging two directions
    [~, rot90] = sw_rot(n,pi/2);
    Sab = (Sab+mmat(mmat(rot90,Sab),rot90'))/2;
    
    % dispersion
    omega = [omega(:,kmIdx==1); omega(:,kmIdx==2); omega(:,kmIdx==3)];
    % exchange matrices
    Sab   = cat(3,mmat(Sab(:,:,:,kmIdx==1),K1), mmat(Sab(:,:,:,kmIdx==2),K2), mmat(Sab(:,:,:,kmIdx==3),conj(K1)));
    
    hkl   = hkl(:,kmIdx==2);
    nHkl0 = nHkl0/3;
end

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
    omega = mat2cell(omega,size(omega,1),repmat(nHkl0,[1 nTwin]));
end

% Creates output structure with the calculated values.
spectra.omega  = omega;
spectra.Sab    = Sab;
spectra.hkl    = hkl(:,1:nHkl0);
spectra.hklA   = hklA;
spectra.incomm = incomm;
% save the important parameters
spectra.param.notwin    = param.notwin;
spectra.param.sortMode  = param.sortMode;
spectra.param.tol       = param.tol;
spectra.param.omega_tol = param.omega_tol;

if ~param.fitmode
    spectra.ff  = ones(nMagExt,nHkl0);
    spectra.obj = copy(obj);
end


end