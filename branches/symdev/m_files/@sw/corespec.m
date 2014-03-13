function [omega, Vsave] = corespec(obj, hkl, varargin)
% calculates dispersion and V transformation matrix using linear spin wave theory
%
% [omega, V] = CORESPEC(obj, k, 'option1', value1 ...)
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
% sortMode      The spin wave modes will be sorted if true. Default is
%               true.
% optmem        Parameter to optimise memory usage. The list of hkl values
%               will be cut into optmem number of pieces and will be
%               calculated piece by piece to decrease memory usage. Default
%               of optmem is zero, when the number of slices are determined
%               automatically from the available free memory.
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
% nSlice        Parameter to optimise memory usage. The list of hkl values
%               will be cut into nSlice pieces and will be calculated piece
%               by piece to decrease memory usage. Default is zero, when
%               the number of slices are determined automatically from the
%               available free memory.
% onlyV         Calculate only the dispersion and the V linear
%               transformation matrices. The V matrices transform between
%               the original magnon operators and the normal magnon
%               operators.
%
% Output:
%
% spectra is a structure, with the following fields:
% omega         Calculated spin wave dispersion, dimensins are
%               [nMode nHkl], where nMagExt is the number of magnetic
%               atoms in the extended unit cell.
% Sab           Dynamical structure factor, dimensins are
%               [3 3 nMode nHkl]. Each (:,:,i,j) submatrix contains the
%               9 correlation functions: Sxx, Sxy, Sxz, etc.
% V             Optional output, calculated if 'onlyV' option is true (in
%               that case Sab is not calculated). Dimensions are
%               [nMode nMode nHkl].
%
%
% nMode is the number of magnetic mode. For commensurate structures it is
% double the number of magnetic atoms in the magnetic cell/supercell. For
% incommensurate structures this number is tripled due to the appearance of
% the (Q+/-km) Fourier components in the correlation functions.
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
% See also SW, SW.SPINWAVESYM, SW_NEUTRON, SW_POL, SW.POWSPEC, SW.OPTMAGSTR.
%

% help when executed without argument
if nargin==1
    help sw.core
    return
end


inpForm.fname  = {'fitmode' 'sortMode' 'optmem' 'fid' 'tol'  'omega_tol' 'hermit' };
inpForm.defval = {false     true       0        0      1e-4  1e-5        true     };
inpForm.size   = {[1 1]     [1 1]      [1 1]    [1 1]  [1 1] [1 1]       [1 1]    };

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

fid = param.fid;

nHkl    = size(hkl,2);

% Create the interaction matrix and atomic positions in the extended
% magnetic unit cell.
if param.fitmode
    [SS, SI, ~] = obj.intmatrix('fitmode',2);
else
    [SS, SI, ~] = obj.intmatrix;
end

% Introduce the opposite couplings.
% (i-->j) and (j-->i)
% transpose the JJ matrix as well [1 2 3 4 5 6 7 8 9] --> [6 9 12 7 10 13 8 11 14]
SS.new         = [SS.all(1:3,:)   -SS.all(1:3,:)  ];
SS.new(4:5,:)  = [SS.all([4 5],:)  SS.all([5 4],:)];
SS.new(6:14,:) = [SS.all(6:14,:)   SS.all([6 9 12 7 10 13 8 11 14],:) ]/2;
SS.all         = SS.new;

% Converts wavevctor list into the extended unit cell
hklExt  = bsxfun(@times,hkl,nExt')*2*pi;

M0 = obj.mag_str.S;
S0 = sqrt(sum(M0.^2,1));

% normal to rotation of the magnetic moments
n  = obj.mag_str.n;
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

% Magnetic field
twinB = SI.field*obj.unit.muB;
MF = repmat(twinB*permute(mmat(SI.g,permute(eta,[1 3 2])),[1 3 2]),[1 2]);

% Creates the serial indices for every matrix element in ham matrix.
idxA1 = [atom1'         atom2'         ];
idxA2 = [atom1'         atom1'         ];
idxB  = [atom1'         atom2'+nMagExt ];
% transpose of idxB
idxC  = [atom2'+nMagExt atom1'         ];
idxD1 = idxA1+nMagExt;
idxD2 = [atom2'+nMagExt atom2'+nMagExt ];
idxMF = [(1:2*nMagExt)' (1:2*nMagExt)' ];

if param.optmem == 0
    freeMem = sw_freemem;
    if freeMem > 0
        nSlice = ceil(nMagExt^2*nHkl*6912/sw_freemem*2);
    else
        nSlice = 1;
        warning('sw:spinwave:FreeMemSize','The size of the free memory is unkown, no memory optimisation!');
    end
else
    nSlice = param.optmem;
end

if nHkl < nSlice
    if fid ~= 0
        fprintf(fid,['Memory allocation is not optimal, nMagExt is'...
            ' too large compared to the free memory!\n']);
    end
    nSlice = nHkl;
elseif nSlice > 1
    if fid ~= 0
        fprintf(fid,['To optimise memory allocation, Q is cut'...
            ' into %d pieces!\n'],nSlice);
    end
end

hklIdx = [floor(((1:nSlice)-1)/nSlice*nHkl)+1 nHkl+1];

% Empty omega dispersion of all spin wave modes, size: 2*nMagExt x nHkl.
omega = zeros(2*nMagExt,0);

if fid == 1
    sw_status(0,1);
end

warn1 = false;

Vsave = zeros(2*nMagExt,2*nMagExt,0);

for jj = 1:nSlice
    % q indices selected for every chunk
    hklIdxMEM  = hklIdx(jj):(hklIdx(jj+1)-1);
    % q values
    hklExtMEM  = hklExt(:,hklIdxMEM);
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
        ham = ham + repmat(accumarray(idxMF,MF,[1 1]*2*nMagExt),[1 1 nHklMEM]);
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
            
            % sort eigenvalues to decreasing order
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
    
    Vsave = cat(3,Vsave,V);
    if fid == 1
        sw_status(jj/nSlice*100);
    end
    
end

if fid == 1
    sw_status(100,2);
else
    if fid ~= 0
        fprintf(fid,'Calculation finished.\n');
    end
end

if warn1
    warning('sw:spinwave:NonPosDefHamiltonian',['To make the Hamiltonian '...
        'positive definite, a small omega_tol value was added to its diagonal!'])
end

end