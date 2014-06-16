function spectra = spinwave(obj, hkl, varargin)
% calculates dynamical spin-spin correlation function using linear spin wave theory
%
% spectra = SPINWAVE(obj, k, 'option1', value1 ...)
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
%               For symbolic calculation at a general reciprocal space
%               point use sym class input. For example to calculate the
%               spectrum along (h,0,0): hkl = [sym('h') 0 0]. To
%               do calculation at a specific point do for example 
%               sym([0 1 0]), to calculate the spectrum at (0,1,0).
%
% Options:
%
% fitmode       Speedup (for fitting mode only), default is false.
% notwin        If true, the spectra of the twins won't be calculated.
%               Default is false.
% sortMode      The spin wave modes will be sorted if true. Default is
%               true.
% optmem        Parameter to optimise memory usage. The list of hkl values
%               will be cut into optmem number of pieces and will be
%               calculated piece by piece to decrease memory usage. Default
%               of optmem is zero, when the number of slices are determined
%               automatically from the available free memory.
% tol           Tolerance of the incommensurability of the magnetic
%               ordering wavevector. Deviations from integer values of the
%               ordering wavevector smaller than the tolerance are
%               considered to be commensurate. Default value is 1e-4.
% omega_tol     Tolerance on the energy difference of degenerate modes when
%               diagonalising the quadratic form, default is 1e-5.
% hermit        Method for matrix diagonalization:
%                   true      J.H.P. Colpa, Physica 93A (1978) 327,
%                   false     R.M. White, PR 139 (1965) A450.
%               Colpa: the grand dynamical matrix is converted into another
%                      Hermitian matrix, that will give the real
%                      eigenvalues.
%               White: the non-Hermitian g*H matrix will be diagonalised,
%                      that is not the elegant method.
%               Advise:
%               Always use Colpa's method, except when small imaginary
%               eigenvalues are expected. In this case only White's method
%               work. The solution in this case is wrong, however by
%               examining the eigenvalues it can give a hint where the
%               problem is.
%               Default is true.
% saveH         If true, the quadratic form of the Hamiltonian is saved. Be
%               carefull, it can take up lots of memory. Default is false.
% saveT         If true, the matrices that transform the normal magnon
%               modes into the magnon modes localized on the spins are
%               saved. Be carefull, it can take up lots of memory.
%               Default is false.
% saveSabp      If true, the dynamical structure factorin the rotating
%               frame is saved S'(k,omega). Default is false.
%
% Output:
%
% 'spectra' is a structure, with the following fields:
% omega         Calculated spin wave dispersion, dimensins are
%               [nMode nHkl], where nMagExt is the number of magnetic
%               atoms in the extended unit cell.
% Sab           Dynamical structure factor, dimensins are
%               [3 3 nMode nHkl]. Each (:,:,i,j) submatrix contains the
%               9 correlation functions: Sxx, Sxy, Sxz, etc.
% H             Quadratic for mof the Hamiltonian.
%               Only saved if saveH is true.
% T             Transformation matrix from the normal magnon modes to the
%               magnons localized on spins:
%                   x_i = sum_j T_ij * x_j'
%               Only saved if saveT is true.
% Sabp          Dynamical structure factor in the rotating frame,
%               dimensions are [3 3 nMode nHkl], but the number of modes
%               are equal to twice the number of magnetic atoms.
%
% nMode is the number of magnetic mode. For commensurate structures it is
% double the number of magnetic atoms in the magnetic cell/supercell. For
% incommensurate structures this number is tripled due to the appearance of
% the (Q+/-km) Fourier components in the correlation functions. For every k
% points in the following order: (k-km,k,k+km).
%
% If several twins exist in the sample, omega and Sab are packaged into a
% cell, that contains nTwin number of matrices.
%
% hkl           Contains the input Q values, dimensins are [3 nHkl].
% hklA          Same Q values, but in reciproc Angstrom units in the
%               lab coordinate system, dimensins are [3 nHkl].
% incomm        Whether the spectra calculated is incommensurate or not.
% obj           The copy of the input obj.
%
% Example:
%
% tri = sw_model('triAF',1);
% sw_plotspec(tri.spinwave({[0 0 0] [1 1 0]}))
%
% The above example will calculate and plot the spin wave dispersion of the
% triangular lattice antiferromagnet (S=1, J=1) along the [H H 0] direction
% in reciprocal space.
%
% See also SW, SW.SPINWAVESYM, SW_NEUTRON, SW.POWSPEC, SW.OPTMAGSTR, SW.FILEID.
%

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

% calculate symbolic spectrum if obj is in symbolic mode
if obj.symb
    if numel(hkl) == 3
        hkl = sym(hkl);
    end
    
    if ~isa(hkl,'sym')
        inpForm.fname  = {'fitmode'};
        inpForm.defval = {false    };
        inpForm.size   = {[1 1]    };
        param0 = sw_readparam(inpForm, varargin{:});
        
        if ~param0.fitmode
            fprintf0(obj.fileid,'No hkl value was given, spin wave spectrum for general Q (h,k,l) will be calculated!');
        end
        spectra = obj.spinwavesym(varargin{:});
    else
        spectra = obj.spinwavesym(varargin{:},'hkl',hkl);
    end
    return
end

% help when executed without argument
if nargin==1
    help sw.spinwave
    spectra = [];
    return
end

inpForm.fname  = {'fitmode' 'notwin' 'sortMode' 'optmem' 'tol' 'hermit'};
inpForm.defval = {false     false    true       0        1e-4  true    };
inpForm.size   = {[1 1]     [1 1]    [1 1]      [1 1]    [1 1] [1 1]   };

inpForm.fname  = [inpForm.fname  {'omega_tol' 'saveSabp' 'saveT' 'saveH'}];
inpForm.defval = [inpForm.defval {1e-5        true       false   false  }];
inpForm.size   = [inpForm.size   {[1 1]       [1 1]      [1 1]   [1 1]  }];

param = sw_readparam(inpForm, varargin{:});

if param.fitmode
    param.sortMode = false;
end

% size of the extended magnetic unit cell
nExt    = double(obj.mag_str.N_ext);
% magnetic ordering wavevector in the extended magnetic unit cell
km = obj.mag_str.k.*nExt;
% whether the structure is incommensurate
incomm = any(abs(km-round(km)) > param.tol);

% Calculates momentum transfer in A^-1 units.
hklA = 2*pi*(hkl'/obj.basisvector)';

% Check for 2*km
tol = param.tol*2;
helical =  sum(abs(mod(abs(2*km)+tol,1)-tol).^2) > tol;

if incomm
   if ~helical
       error('sw:spinwave:Twokm',['The two times the magnetic ordering '...
           'wavevector 2*km = G, reciproc lattice vector, use magnetic supercell to calculate spectrum!']);
   end
    % without the k_m: (k, k, k)
    hkl0 = repmat(hkl,[1 3]);
    % calculate dispersion for (k-km, k, k+km)
    hkl  = [bsxfun(@minus,hkl,km') hkl bsxfun(@plus,hkl,km')];
else
    hkl0 = hkl;
end

% Print output into the following file
fid = obj.fid;

spectra = struct;
nHkl    = size(hkl,2);


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
    % without the +/- k_m value
    hkl0 = cell2mat(obj.twinq(hkl0));
    nHkl = nHkl*nTwin;
end

% determines a twin index for every q point
twinIdx = repmat(1:nTwin,[nHkl 1]);
twinIdx = twinIdx(:);

% Create the interaction matrix and atomic positions in the extended
% magnetic unit cell.
if param.fitmode
    [SS, SI, RR] = obj.intmatrix('fitmode',2,'conjugate',true);
else
    [SS, SI, RR] = obj.intmatrix('conjugate',true);
end

% Converts wavevctor list into the extended unit cell
hklExt  = bsxfun(@times,hkl,nExt')*2*pi;
% q values without the +/-k_m value
hklExt0 = bsxfun(@times,hkl0,nExt')*2*pi;

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
        fprintf0(fid,['Calculating INCOMMENSURATE spin wave spectra '...
            '(nMagExt = %d, nHkl = %d, nTwin = %d)...\n'],nMagExt, nHkl0, nTwin);
    else
        fprintf0(fid,['Calculating COMMENSURATE spin wave spectra '...
            '(nMagExt = %d, nHkl = %d, nTwin = %d)...\n'],nMagExt, nHkl0, nTwin);
    end
end

% Local (e1,e2,e3) coordinate system fixed to the moments, 
% e3||Si, 
% e2 = Si x [1,0,0], if Si || [1,0,0] --> e2 = [0,0,1]
% e1 = e2 x e3
magTab = obj.magtable;

zed = magTab.e1 + 1i*magTab.e2;
eta = magTab.e3;

dR    = [SS.all(1:3,:) zeros(3,nMagExt)];
atom1 = [SS.all(4,:)   1:nMagExt];
atom2 = [SS.all(5,:)   1:nMagExt];
% magnetic couplings, 3x3xnJ
JJ = cat(3,reshape(SS.all(6:end,:),3,3,[]),SI.aniso);

if incomm
    % transform JJ due to the incommensurate wavevector
    [~, K] = sw_rot(n,km*dR*2*pi);
    % multiply JJ with K matrices for every interaction
    % and symmetrising JJ for the rotating basis
    JJ = (mmat(JJ,K)+mmat(K,JJ))/2;
end

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

% Magnetic field is different for every twin
%MF  =  repmat(obj.unit.muB*SI.field*eta,[1 2]);
MF = zeros(1,2*nMagExt,nTwin);
for ii = 1:nTwin
    % rotate the magnetic field to the relative direction of every twin
    % backward rotation with the rotc matrix of the twin
    twinB = SI.field*obj.twin.rotc(:,:,ii)*obj.unit.muB;
    MF(:,:,ii) = repmat(twinB*permute(mmat(SI.g,permute(eta,[1 3 2])),[1 3 2]),[1 2]);
end

% Creates the serial indices for every matrix element in ham matrix.
idxA1 = [atom1'         atom2'         ];
idxA2 = [atom1'         atom1'         ];
idxB  = [atom1'         atom2'+nMagExt ];
% transpose of idxB
idxC  = [atom2'+nMagExt atom1'         ];
idxD1 = idxA1+nMagExt;
idxD2 = [atom2'+nMagExt atom2'+nMagExt ];
idxMF = [(1:2*nMagExt)' (1:2*nMagExt)' ];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MEMORY MANAGEMENT LOOP
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if param.optmem == 0
    freeMem = sw_freemem;
    if freeMem > 0
        nSlice = ceil(nMagExt^2*nHkl*6912/sw_freemem*2);
    else
        nSlice = 1;
        if ~param.fitmode
            warning('sw:spinwave:FreeMemSize','The size of the free memory is unkown, no memory optimisation!');
        end
    end
else
    nSlice = param.optmem;
end

if nHkl < nSlice
    if fid ~= 0
        fprintf0(fid,['Memory allocation is not optimal, nMagExt is'...
            ' too large compared to the free memory!\n']);
    end
    nSlice = nHkl;
elseif nSlice > 1
    if fid ~= 0
        fprintf0(fid,['To optimise memory allocation, Q is cut'...
            ' into %d pieces!\n'],nSlice);
    end
end

hklIdx = [floor(((1:nSlice)-1)/nSlice*nHkl)+1 nHkl+1];

% Empty omega dispersion of all spin wave modes, size: 2*nMagExt x nHkl.
omega = zeros(2*nMagExt,0);

% empty Sab
Sab = zeros(3,3,2*nMagExt,0);

% Empty matrices to save different intermediate results for further
% analysis: Hamiltonian, eigenvectors, dynamical structure factor in the
% rotating frame
if param.saveT
    Tsave = zeros(2*nMagExt,2*nMagExt,nHkl);
end
if param.saveH
    Hsave = zeros(2*nMagExt,2*nMagExt,nHkl);
end

if fid == 1
    sw_status(0,1);
end

warn1 = false;

for jj = 1:nSlice
    % q indices selected for every chunk
    hklIdxMEM  = hklIdx(jj):(hklIdx(jj+1)-1);
    % q values contatining the k_m vector
    hklExtMEM  = hklExt(:,hklIdxMEM);
    % q values without the +/-k_m vector
    hklExt0MEM = hklExt0(:,hklIdxMEM);
    % twin indices for every q point
    twinIdxMEM = twinIdx(hklIdxMEM);
    nHklMEM = size(hklExtMEM,2);
    
    % Creates the matrix of exponential factors nCoupling x nHkl size.
    % Extends dR into 3 x 3 x nCoupling x nHkl
    ExpF = exp(1i*permute(sum(repmat(dR,[1 1 nHklMEM]).*repmat(...
        permute(hklExtMEM,[1 3 2]),[1 nCoupling 1]),1),[2 3 1]))';
    
    % Creates the matrix elements containing zed.
    A1 = bsxfun(@times,     AD0 ,ExpF);
    B  = bsxfun(@times,     BC0 ,ExpF);
    D1 = bsxfun(@times,conj(AD0),ExpF);
    
    % Store all indices
    idxAll = [idxA1; idxB; idxC; idxD1];
    % Store all matrix elements
    ABCD   = [A1     B     conj(B)  D1];
    
    % Stores the matrix elements in ham.
    idx3   = repmat(1:nHklMEM,[4*nCoupling 1]);
    idxAll = [repmat(idxAll,[nHklMEM 1]) idx3(:)];
    idxAll = idxAll(:,[2 1 3]);
    
    ABCD   = ABCD';
    
    ham = accumarray(idxAll,ABCD(:),[2*nMagExt 2*nMagExt nHklMEM]);
    
    ham = ham + repmat(accumarray([idxA2; idxD2],2*[A20 D20],[1 1]*2*nMagExt),[1 1 nHklMEM]);
    
    if any(SI.field)
        % different field for different twin
        for ii = min(twinIdxMEM):max(twinIdxMEM)
            nTwinQ = sum(twinIdxMEM==ii);
            ham(:,:,twinIdxMEM==ii) = ham(:,:,twinIdxMEM==ii) + ...
                repmat(accumarray(idxMF,MF(:,:,ii),[1 1]*2*nMagExt),[1 1 nTwinQ]);
        end
        
        %ham = ham + repmat(accumarray(idxMF,MF,[1 1]*2*nMagExt),[1 1 nHklMEM]);
    end
    
    ham = (ham + conj(permute(ham,[2 1 3])))/2;
    
    g  = diag([ones(nMagExt,1); -ones(nMagExt,1)]);
    gd = diag(g);
    
    if param.hermit
        % All the matrix calculations are according to Colpa's paper
        % J.H.P. Colpa, Physica 93A (1978) 327-353
        
        V = zeros(2*nMagExt,2*nMagExt,nHklMEM);
        
        for ii = 1:nHklMEM
            [K, posDef]  = chol(ham(:,:,ii));
            if posDef > 0
                try
                    K = chol(ham(:,:,ii)+eye(2*nMagExt)*param.omega_tol);
                    warn1 = true;
                catch PD
                    error('sw:spinwave:NonPosDefHamiltonian',...
                        ['Hamiltonian matrix is not positive definite, probably'...
                        ' the magnetic structure is wrong! For approximate'...
                        ' diagonalization try the param.hermit=false option']);
                end
            end
            
            K2 = K*g*K';
            K2 = 1/2*(K2+K2');
            % Hermitian K2 will give orthogonal eigenvectors
            if param.sortMode
                [U, D] = eigenshuffle(K2);
            else
                [U, D] = eig(K2);
                D = diag(D);
            end
            
            % sort eigenvalues to decreasing order this contradicts with
            % eigenshuffle
            
            % TODO
            [D, idx] = sort(D,'descend');
            U = U(:,idx);
            
            % omega dispersion
            omega(:,end+1) = D; %#ok<AGROW>
            
            % the inverse of the para-unitary transformation V
            V(:,:,ii) = inv(K)*U*diag(sqrt(gd.*omega(:,end))); %#ok<MINV>
        end
    else
        % All the matrix calculations are according to White's paper
        % R.M. White, et al., Physical Review 139, A450?A454 (1965)
        
        gham = 0*ham;
        for ii = 1:nHklMEM
            gham(:,:,ii) = g*ham(:,:,ii);
        end
        
        [V, D] = eigorth(gham,param.omega_tol, param.sortMode);
        
        for ii = 1:nHklMEM
            % multiplication with g removed to get negative and positive
            % energies as well
            omega(:,end+1) = D(:,ii); %#ok<AGROW>
            M           = diag(g*V(:,:,ii)'*g*V(:,:,ii));
            V(:,:,ii)   = V(:,:,ii)*diag(sqrt(1./M));
        end
    end
    
    if param.saveT
        Tsave(:,:,hklIdxMEM) = V;
    end
    if param.saveH
        Hsave(:,:,hklIdxMEM) = ham;
    end
    
    % Calculates correlation functions.
    % V right
    VExtR = repmat(permute(V  ,[4 5 1 2 3]),[3 3 1 1 1]);
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
    
    if fid == 1
        sw_status(jj/nSlice*100);
    end
    
end

if fid == 1
    sw_status(100,2);
else
    if fid ~= 0
        fprintf0(fid,'Calculation finished.\n');
    end
end

if warn1 && ~param.fitmode
    warning('sw:spinwave:NonPosDefHamiltonian',['To make the Hamiltonian '...
        'positive definite, a small omega_tol value was added to its diagonal!'])
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
    
    % keep the rotation invariant part of Sab
    %nx  = [0 -n(3) n(2);n(3) 0 -n(1);-n(2) n(1) 0];
    %nxn = n'*n;
    m1  = eye(3);
    
    % if the 2*km vector is integer, the magnetic structure is not a true
    % helix
    %tol = param.tol*2;
    %helical =  sum(abs(mod(abs(2*km)+tol,1)-tol).^2) > tol;
    
    if helical
        % integrating out the arbitrary initial phase of the helix
        Sab = 1/2*Sab - 1/2*mmat(mmat(nx,Sab),nx) + 1/2*mmat(mmat(nxn-m1,Sab),nxn) + 1/2*mmat(mmat(nxn,Sab),2*nxn-m1);
    end
    
    % Save the structure factor in the rotating frame
    if param.saveSabp
        Sabp = Sab(:,:,:,kmIdx==2);
        omegap = omega(:,kmIdx==2);
    end
    
    % dispersion
    omega = [omega(:,kmIdx==1); omega(:,kmIdx==2); omega(:,kmIdx==3)];
    % exchange matrices
    Sab   = cat(3,mmat(Sab(:,:,:,kmIdx==1),K1), mmat(Sab(:,:,:,kmIdx==2),K2), ...
        mmat(Sab(:,:,:,kmIdx==3),conj(K1)));
    
    hkl   = hkl(:,kmIdx==2);
    nHkl0 = nHkl0/3;
else
    helical = false;
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
spectra.omega   = omega;
spectra.Sab     = Sab;
spectra.hkl     = hkl(:,1:nHkl0);
spectra.hklA    = hklA;
spectra.incomm  = incomm;
spectra.helical = helical;
spectra.norm    = false;

% Save different intermediate results.
if param.saveT
    spectra.T = Tsave;
end
if param.saveH
    spectra.H = Hsave;
end
if param.saveSabp && incomm
    spectra.Sabp = Sabp;
    spectra.omegap = omegap;
end

% save the important parameters
spectra.param.notwin    = param.notwin;
spectra.param.sortMode  = param.sortMode;
spectra.param.tol       = param.tol;
spectra.param.omega_tol = param.omega_tol;
spectra.param.hermit    = param.hermit;

if ~param.fitmode
    spectra.ff  = ones(nMagExt,nHkl0);
    spectra.obj = copy(obj);
end

end