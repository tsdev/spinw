function spectra = spinwavesym(obj, varargin)
% calculates symbolic spin wave dispersion
%
% spectra = SPINWAVESYM(obj, 'option1', value1 ...)
%
% Symbolic spin wave dispersion  is calculated as a function of reciprocal
% space points. The function can deal with arbitrary magnetic structure and
% magnetic interactions as well as single ion anisotropy and magnetic
% field.
%
% If the magnetic ordering wavevector is non-integer, the dispersion is
% calculated using a coordinate system rotating from cell to cell. In this
% case the Hamiltonian has to fulfill this extra rotational symmetry.
%
% The method for matrix diagonalization is according to R.M. White, PR 139
% (1965) A450. The non-Hermitian g*H matrix will be diagonalised.
%
% Input:
%
% obj           Input structure, sw class object.
%
% Options:
%
% hkl           Symbolic definition of q vector. Default is the general Q
%               point:
%                   hkl = [sym('h') sym('k') sym('l')]
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

hkl0 = [sym('h','real'); sym('k','real'); sym('l','real')];

inpForm.fname  = {'fitmode' 'notwin' 'sortMode' 'optmem' 'fid' 'tol' 'omega_tol' 'hkl'};
inpForm.defval = {false     false    true      true     1     1e-4   1e-5        hkl0  };
inpForm.size   = {[1 1]     [1 1]    [1 1]      [1 1]    [1 1] [1 1]  [1 1]      [3 1]};

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
%syms h k l real
%hkl = [h; k; l];
hkl = param.hkl;

fid = param.fid;


% Create the interaction matrix and atomic positions in the extended
% magnetic unit cell.
%[SS, SI, RR] = obj.intmatrix('plotmode',true,'extend',true,'fitmode',2);
[SS, SI] = obj.intmatrix('plotmode',true,'extend',true,'fitmode',2,'conjugate',true);

% Introduce the opposite couplings.
% (i-->j) and (j-->i)
% transpose the JJ matrix as well [1 2 3 4 5 6 7 8 9] --> [6 9 12 7 10 13 8 11 14]
% SS.new         = [SS.all(1:3,:)   -SS.all(1:3,:)  ];
% SS.new(4:5,:)  = [SS.all([4 5],:)  SS.all([5 4],:)];
% SS.new(6:14,:) = [SS.all(6:14,:)   SS.all([6 9 12 7 10 13 8 11 14],:) ]/2;

% ???
%idxM = [SS.all(15,:) SS.all(15,:)];

%SS.all         = SS.new;

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
e3 = M0./repmat(S0(:),[3 1]);
% e2 = Si x [1,0,0], if Si || [1,0,0] --> e2 = [0,0,1]
e2  = [zeros(1,nMagExt); e3(3,:); -e3(2,:)];
%e2(3,~any(e2)) = 1;
e2(3,~any(abs(e2)>1e-10)) = 1;
E0 = sqrt(sum(e2.^2,1));
e2  = e2./repmat(E0(:),[3 1]);
% e1 = e2 x e3
e1  = cross(e2,e3);

% Defines eta and zed.
zed = e1 + 1i*e2;
eta = e3;

dR    = [SS.all(1:3,:) zeros(3,nMagExt)];
atom1 = int32([SS.all(4,:)   1:nMagExt]);
atom2 = int32([SS.all(5,:)   1:nMagExt]);
%idxM  = int32([idxM repmat(obj.single_ion.aniso,1,prod(nMagExt))]);
% magnetic couplings, 3x3xnJ
JJ = cat(3,reshape(SS.all(6:14,:),3,3,[]),SI.aniso);

% remove zero anisotropy matrices
anyIdx = squeeze(sumsym(sumsym(abs(JJ),1),2)) == 0;
if ~isa(anyIdx,'logical')
    anyIdx = isAlways(anyIdx);
end
anyIdx = ~anyIdx;

dR    = dR(:,anyIdx);
atom1 = atom1(1,anyIdx);
atom2 = atom2(1,anyIdx);
JJ    = JJ(:,:,anyIdx);
%idxM  = idxM(1,anyIdx);

% create symbolic variables from matrix labels
%nMat = size(obj.matrix.mat,3);

%for ii = 1:nMat
%    JS(ii) = sym(obj.matrix.label{ii},'real'); %#ok<AGROW>
%end

% normalize JJ to the maximum value
%JJ = JJ./repmat(max(max(abs(JJ))),[3 3 1]);


% create symbolic matrices
%JJ = JJ.*repmat(permute(JS(idxM),[3 1 2]),[3 3 1]);

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
A0 =  SiSj.*shiftdim(sumsym(sumsym(zedL.*JJ.*conj(zedR),2),1),1);
B0 =  SiSj.*shiftdim(sumsym(sumsym(zedL.*JJ.*     zedR ,2),1),1);
C0 =  shiftdim(sumsym(sumsym(etaL.*JJ.*etaR,2),1),1);
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

% All the matrix calculations are according to White's paper
% R.M. White, et al., Physical Review 139, A450?A454 (1965)

fprintf(fid,'Calculating SYMBOLIC eigenvalues... ');
[V, D] = eig(g*ham); % 3rd output P
fprintf('ready!\n');

spectra.V0 = V;

if size(V,2) == size(ham,2)
    M = diag(g*V'*g*V);
    V = V*diag(sqrt(1./M));
    spectra.V = V;
else
    warning('There are degenerate eigenvalues!')
end

% multiplication with g removed to get negative and positive
% energies as well
spectra.omega = simplify(diag(D));

end