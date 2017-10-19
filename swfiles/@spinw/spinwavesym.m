function spectra = spinwavesym(obj, varargin)
% calculates symbolic spin wave dispersion
% 
% ### Syntax
% 
% `spectra = spinwavesym(obj,Name,Value)`
% 
% ### Description
% 
% `spectra = spinwavesym(obj,Name,Value)` calculates symbolic spin wave
% dispersion as a function of $Q$. The function can deal with arbitrary
% magnetic structure and magnetic interactions as well as single ion
% anisotropy and magnetic field. Biquadratic exchange interactions are also
% implemented, however only for $k=0$ magnetic structures.
% 
% If the magnetic propagation vector is non-integer, the dispersion is
% calculated using a coordinate system rotating from cell to cell. In this
% case the Hamiltonian has to fulfill this extra rotational symmetry.
%  
% The method works for incommensurate structures, however the calculated
% omega dispersion does not contain the $\omega(\mathbf{k}\pm \mathbf{k}_m)$ terms that has to be
% added manually.
%  
% The method for matrix diagonalization is according to R.M. White, PR 139
% (1965) A450. The non-Hermitian g*H matrix will be diagonalised.
%  
% ### Examples
%
% The first section of the example calculates the symbolic spin wave
% spectrum. Unfortunatelly the symbolic expression needs manipulations to
% bring it to readable form. To check the solution, the second section
% converts the symbolic expression into a numerical vector and the third
% section plots the real and imaginary part of the solution.
%
% ```
% >>tri = sw_model('triAF',1)
% >>tri.symbolic(true)
% >>tri.genmagstr('mode','direct','k',[1/3 1/3 0],'S',[1 0 0])
% >>symSpec = tri.spinwave
% >>pretty(symSpec.omega)>>
% >>J_1 = 1
% >>h = linspace(0,1,500)
% >>k = h
% >>omega = eval(symSpec.omega)
% >>p1 = plot(h,real(omega(1,:)),'k-')
% >>hold on
% >>plot(h,real(omega(2,:)),'k-')
% >>p2 = plot(h,imag(omega(1,:)),'r-')
% >>plot(h,imag(omega(2,:)),'r-')
% >>xlabel('Momentum (h,h,0) (r.l.u.)')
% >>ylabel('Energy (meV)')
% >>legend([p1 p2],'Real(\omega(Q))','Imag(\omega(Q))')
% >>title('Spin wave dispersion of the TLHAF')
% >>snapnow
% ```
%
% ### Input Arguments
% 
% `obj`
% : [spinw] object.
% 
% ### Name-Value Pair Arguments
% 
% `'hkl'`
% : Symbolic definition of $Q$ vector. Default is the general $Q$
%   point:
%   ```
%   hkl = [sym('h') sym('k') sym('l')]
%   ```
% 
% `'eig'`
% : If true the symbolic Hamiltonian is diagonalised symbolically. For
%   large matrices (many magnetic atom per unit cell) this might be
%   impossible. Set `eig` to `false` to output only the quadratic
%   Hamiltonian. Default is `true`.
% 
% `'vect'`
% : If `true` the eigenvectors are also calculated. Default is `false`.
% 
% `'tol'`
% : Tolerance of the incommensurability of the magnetic
%   ordering wavevector. Deviations from integer values of the
%   ordering wavevector smaller than the tolerance are
%   considered to be commensurate. Default value is $10^{-4}$.
% 
% `'norm'`
% : Whether to produce the normalized symbolic eigenvectors. It can be
%   impossible for large matrices, in that case set it to
%   `false`. Default is `true`.
% 
% `'fid'`
% : Defines whether to provide text output. The default value is determined
%   by the `fid` preference stored in [swpref]. The possible values are:
%   * `0`   No text output is generated.
%   * `1`   Text output in the MATLAB Command Window.
%   * `fid` File ID provided by the `fopen` command, the output is written
%           into the opened file stream.
%
% `'title'`
% : Gives a title string to the simulation that is saved in the
%   output.
%
% ### Output Arguments
%
% `spectra`
% : Structure, with the following fields:
%   * `omega`   Calculated spin wave dispersion, dimensins are
%               $[2*n_{magExt}\times n_{hkl}]$, where $n_{magExt}$ is the number of magnetic
%               atoms in the extended unit cell.
%   * `V0`      Eigenvectors of the quadratic Hamiltonian.
%   * `V`       Normalized eigenvectors of the quadratic Hamiltonian.
%   * `ham`     Symbolic matrix of the Hamiltonian.
%   * `incomm`  Whether the spectra calculated is incommensurate or not.
%   * `obj`     The clone of the input `obj`.
%
% ### See Also
 %
% [spinw] \| [spinw.spinwave]
%

% save the begining time of the calculation
spectra.datestart = datestr(now);

hkl0 = [sym('h','real'); sym('k','real'); sym('l','real')];

title0 = 'Symbolic LSWT spectrum';

inpForm.fname  = {'tol' 'hkl'  'eig' 'vect' 'norm' 'title' 'fid'};
inpForm.defval = {1e-4   hkl0   true false  false   title0 -1   };
inpForm.size   = {[1 1] [3 1]  [1 1] [1 1]  [1 1]  [1 -1]  [1 1]};

param = sw_readparam(inpForm, varargin{:});

if param.norm
    param.vect = true;
end

% generate the magnetic structure
magstr = obj.magstr;

% size of the extended magnetic unit cell
nExt    = magstr.N_ext;
% magnetic ordering wavevector in the extended magnetic unit cell
km = magstr.k.*nExt;
% whether the structure is incommensurate
incomm = any(~sw_always(abs(km-round(km)) <= param.tol));

% symbolic wavevectors, convert to the model rlu
hkl = obj.unit.qmat*param.hkl;

if param.fid == -1
    fid = swpref.getpref('fid',true);
else
    fid = param.fid;
end

% Create the interaction matrix and atomic positions in the extended
% magnetic unit cell.
%[SS, SI, RR] = obj.intmatrix('plotmode',true,'extend',true,'fitmode',2);
%if obj.symmetry && any(sw_mattype(obj.matrix.mat)~=1)
%    warning('spinw:spinwavesym:symmetry',['The non-isotropic symbolic matrices will'...
%    'not be rotated unsing the point group operators, define them manually!'])
%end
%[SS, SI] = obj.intmatrix('plotmode',true,'extend',true,'fitmode',2,'conjugate',true,'rotMat',false);
%[SS, SI] = obj.intmatrix('plotmode',true,'extend',true,'fitmode',2,'conjugate',true);
[SS, SI] = obj.intmatrix('fitmode',true,'conjugate',true);

% is there any biquadratic exchange
bq = SS.all(15,:)==1;

% Biquadratic exchange only supported for commensurate structures
if incomm && any(bq)
    error('spinw:spinwavesym:Biquadratic','Biquadratic exchange can be only calculated for k=0 structures!');
end

if any(bq)
    % Separate the biquadratic couplings
    % Just use the SS.bq matrix produced by intmatrix(), it won't contain
    % the transpose matrices (not necessary for biquadratic exchange)
    % TODO check whether to keep the transposed matrices to be sure
    SS.bq = SS.all(1:6,bq);
    % Keep only the quadratic exchange couplings
    SS.all = SS.all(1:14,SS.all(15,:)==0);
end

% Converts wavevctor list into the extended unit cell
hklExt  = hkl.*nExt'*2*pi;

% Calculates parameters eta and zed.
if isempty(magstr.S)
    error('spinw:spinwave:NoMagneticStr','No magnetic structure defined in obj!');
end

M0 = magstr.S;
S0 = sqrt(sum(M0.^2,1));
% normal to rotation of the magnetic moments
n  = magstr.n;
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
anyIdx = ~sw_always(squeeze(sumsym(sumsym(abs(JJ),1),2)) == 0);

dR    = dR(:,anyIdx);
atom1 = atom1(1,anyIdx);
atom2 = atom2(1,anyIdx);
JJ    = JJ(:,:,anyIdx);


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

% Calculate matrix elements for biquadratic exchange
if any(bq)
    bqdR    = SS.bq(1:3,:);
    bqAtom1 = int32(SS.bq(4,:));
    bqAtom2 = int32(SS.bq(5,:));
    bqJJ    = SS.bq(6,:);

    % matrix elements: M,N,P,Q
    bqM = sum(eta(:,bqAtom1).*eta(:,bqAtom2),1);
    bqN = sum(eta(:,bqAtom1).*zed(:,bqAtom2),1);
    bqO = sum(zed(:,bqAtom1).*zed(:,bqAtom2),1);
    bqP = sum(conj(zed(:,bqAtom1)).*zed(:,bqAtom2),1);
    bqQ = sum(zed(:,bqAtom1).*eta(:,bqAtom2),1);
    
    Si = S0(bqAtom1);
    Sj = S0(bqAtom2);
    % C_ij matrix elements
    bqA0 = (Si.*Sj).^(3/2).*(bqM.*conj(bqP) + bqQ.*conj(bqN)).*bqJJ;
    bqB0 = (Si.*Sj).^(3/2).*(bqM.*bqO + bqQ.*bqN).*bqJJ;
    bqC  = Si.*Sj.^2.*(conj(bqQ).*bqQ - 2*bqM.^2).*bqJJ;
    bqD  = Si.*Sj.^2.*(bqQ).^2.*bqJJ;

    % Creates the serial indices for every matrix element in ham matrix.
    % Aij(k) matrix elements (b^+ b)
    idxbqA  = [bqAtom1' bqAtom2'];
    % b b^+ elements
    idxbqA2 = [bqAtom1' bqAtom2']+nMagExt;
    
    % Bij(k) matrix elements (b^+ b^+)
    idxbqB  = [bqAtom1' bqAtom2'+nMagExt];
    % transpose of B (b b)
    %idxbqB2 = [bqAtom2'+nMagExt bqAtom1']; % SP2
    
    idxbqC  = [bqAtom1' bqAtom1'];
    idxbqC2 = [bqAtom1' bqAtom1']+nMagExt;
    
    idxbqD  = [bqAtom1' bqAtom1'+nMagExt];
    %idxbqD2 = [bqAtom1'+nMagExt bqAtom1]; % SP2

end

ham = sym(zeros(2*nMagExt));

for ii = 1:size(idxAll,1)
    ham(idxAll(ii,1),idxAll(ii,2)) = ham(idxAll(ii,1),idxAll(ii,2)) + ABCD(ii);
end

if any(bq)
    bqExpF = exp(1i*hklExt'*bqdR);

    bqA  = bqA0.*bqExpF;
    bqA2 = conj(bqA0).*bqExpF;
    bqB  = bqB0.*bqExpF;

    idxbqAll = [idxbqA; idxbqA2; idxbqB; idxbqC; idxbqC2; idxbqD];
    bqABCD = [bqA bqA2 2*bqB bqC bqC 2*bqD];
    
    for ii = 1:size(idxbqAll,1)
        ham(idxbqAll(ii,1),idxbqAll(ii,2)) = ham(idxbqAll(ii,1),idxbqAll(ii,2)) + bqABCD(ii);
    end
end

ham = simplify((ham + conj(permute(ham,[2 1 3])))/2);

% save symbolic Hamiltonian matrix
spectra.ham = ham;

if param.eig
    g = sym(diag([ones(nMagExt,1); -ones(nMagExt,1)]));
    
    % All the matrix calculations are according to White's paper
    % R.M. White, et al., Physical Review 139, A450?A454 (1965)
    
    fprintf0(fid,'Calculating SYMBOLIC eigenvalues... ');
    
    if param.vect
        [V, D] = eig(g*ham); % 3rd output P
        D = diag(D);
    else
        D = eig(g*ham);
    end
    
    fprintf0(fid,'ready!\n');
    
    if param.vect
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
    end
    
    % multiplication with g removed to get negative and positive
    % energies as well
    spectra.omega = simplify(D);
end

spectra.obj      = copy(obj);
spectra.dateend  = datestr(now);
spectra.title    = param.title;

end