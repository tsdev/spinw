function spectra = spinwavesym(obj, varargin)
% calculates symbolic spin wave dispersion
%
% spectra = SPINWAVESYM(obj, 'option1', value1 ...)
%
% Symbolic spin wave dispersion is calculated as a function of reciprocal
% space points. The function can deal with arbitrary magnetic structure and
% magnetic interactions as well as single ion anisotropy and magnetic
% field.
%
% If the magnetic ordering wavevector is non-integer, the dispersion is
% calculated using a coordinate system rotating from cell to cell. In this
% case the Hamiltonian has to fulfill this extra rotational symmetry.
%
% The method works for incommensurate structures, however the calculated
% omega dispersion does not contain the omega(k+/-km) terms that has to be
% added manually.
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
% eig           If true the symbolic Hamiltonian is diagonalised. For large
%               matrices (many magnetic atom per unit cell) this might be
%               impossible. Set 'eig' to false to output only the quadratic
%               Hamiltonian. Default is true.
% tol           Tolerance of the incommensurability of the magnetic
%               ordering wavevector. Deviations from integer values of the
%               ordering wavevector smaller than the tolerance are
%               considered to be commensurate. Default value is 1e-4.
% norm          Whether to produce the normalized eigenvectors. It can be
%               impossible for large matrices, in that case set it to
%               false. Default is true.
% title         Gives a title string to the simulation that is saved in the
%               output.

% Output:
%
% 'spectra' is a structure, with the following fields:
% omega         Calculated spin wave dispersion, dimensins are
%               [2*nMagExt nHkl], where nMagExt is the number of magnetic
%               atoms in the extended unit cell.
% V0            Eigenvectors of the quadratic Hamiltonian.
% V             Normalized eigenvectors of the quadratic Hamiltonian.
% ham           Symbolic matrix of the Hamiltonian.
%
% If several domains exist in the sample, omega and Sab are packaged into a
% cell, that contains nTwin number of matrices.
%
% incomm        Whether the spectra calculated is incommensurate or not.
% obj           The copy of the input obj.
%
% Example:
%
% tri = sw_model('triAF',1);
% tri.symbolic(true)
% tri.genmagstr('mode','direct','k',[1/3 1/3 0],'S',[1 0 0])
% symSpec = tri.spinwave;
%
% J1 = 1;
% h = linspace(0,1,500);
% k = h;
% omega = eval(symSpec.omega);
%
% p1 = plot(h,real(omega(1,:)),'ko');
% hold on
% plot(h,real(omega(2,:)),'ko')
% p2 = plot(h,imag(omega(1,:)),'r-');
% plot(h,imag(omega(2,:)),'r-')
% xlabel('Momentum (H,H,0) (r.l.u.)')
% ylabel('Energy (meV)')
% legend([p1 p2],'Real(\omega(Q))','Imag(\omega(Q))')
% title('Spin wave dispersion of the TLHAF')
%
% The first section calculates the symbolic spin wave spectrum.
% Unfortunatelly the symbolic expression needs manipulations to bring it to
% readable form. To check the solution, the second section converts the
% symbolic expression into a numerical vector and the third section plots
% the real and imaginary part of the solution.
%
% See also SW, SW.SPINWAVE, SW_NEUTRON, SW_POL, SW.POWSPEC, SW.OPTMAGSTR.
%

% save the begining time of the calculation
spectra.datestart = datetime;

hkl0 = [sym('h','real'); sym('k','real'); sym('l','real')];

title0 = 'Symbolical LSWT spectrum';

inpForm.fname  = {'tol' 'hkl'  'eig' 'norm' 'title'};
inpForm.defval = {1e-4   hkl0   true true   title0 };
inpForm.size   = {[1 1] [3 1]  [1 1] [1 1]  [1 -1] };

param = sw_readparam(inpForm, varargin{:});

% seize of the extended magnetic unit cell
nExt    = double(obj.mag_str.N_ext);
% magnetic ordering wavevector in the extended magnetic unit cell
km = obj.mag_str.k.*nExt;
% whether the structure is incommensurate
incomm = abs(km-round(km)) <= param.tol;
if islogical(incomm)
    incomm = any(~incomm);
else
    incomm = any(~isAlways(incomm));
end

% symbolic wavevectors
hkl = param.hkl;

fid = obj.fid;


% Create the interaction matrix and atomic positions in the extended
% magnetic unit cell.
%[SS, SI, RR] = obj.intmatrix('plotmode',true,'extend',true,'fitmode',2);
%if obj.symmetry && any(sw_mattype(obj.matrix.mat)~=1)
%    warning('sw:spinwavesym:symmetry','The non-isotropic symbolic matrices will not be rotated unsing the point group operators, define them manually!')
%end
%[SS, SI] = obj.intmatrix('plotmode',true,'extend',true,'fitmode',2,'conjugate',true,'rotMat',false);
[SS, SI] = obj.intmatrix('plotmode',true,'extend',true,'fitmode',2,'conjugate',true);

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
        fprintf0(fid,'Calculating SYMBOLIC INCOMMENSURATE spin wave spectra (nMagExt = %d)...\n',nMagExt);
    else
        fprintf0(fid,'Calculating SYMBOLIC COMMENSURATE spin wave spectra (nMagExt = %d)...\n',nMagExt);
    end
end

% Local (e1,e2,e3) coordinate system fixed to the moments,
% e3||Si,
% e2 = Si x [1,0,0], if Si || [1,0,0] --> e2 = [0,0,1]
% e1 = e2 x e3
magTab = obj.magtable;

zed = magTab.e1 + 1i*magTab.e2;
eta = magTab.e3;

if numel(SS.all) > 0
    dR    = [SS.all(1:3,:) zeros(3,nMagExt)];
    atom1 = int32([SS.all(4,:)   1:nMagExt]);
    atom2 = int32([SS.all(5,:)   1:nMagExt]);
    %idxM  = int32([idxM repmat(obj.single_ion.aniso,1,prod(nMagExt))]);
    % magnetic couplings, 3x3xnJ
    JJ = cat(3,reshape(SS.all(6:14,:),3,3,[]),SI.aniso);
else
    % no magnetic couplings
    dR    = zeros(3,nMagExt);
    atom1 = int32(1:nMagExt);
    atom2 = int32(1:nMagExt);
    
    JJ = SI.aniso;
end

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

% save symbolic Hamiltonian matrix
spectra.ham = ham;

if param.eig
    g = sym(diag([ones(nMagExt,1); -ones(nMagExt,1)]));
    
    % All the matrix calculations are according to White's paper
    % R.M. White, et al., Physical Review 139, A450?A454 (1965)
    
    fprintf0(fid,'Calculating SYMBOLIC eigenvalues... ');
    [V, D] = eig(g*ham); % 3rd output P
    fprintf0(fid,'ready!\n');
    
    
    spectra.V0 = V;
    
    if size(V,2) == size(ham,2)
        if param.norm
            % normalized iegenvectors
            M = diag(g*V'*g*V);
            V = V*diag(sqrt(1./M));
            spectra.V = V;
        end
    else
        warning('There are degenerate eigenvalues!')
    end
    
    % multiplication with g removed to get negative and positive
    % energies as well
    spectra.omega = simplify(diag(D));
end

spectra.obj      = copy(obj);
spectra.dateend  = datenum;
spectra.title    = param.title;

end