function M = moment(obj, varargin)
% calculates the size of the reduced moment due to quantum and thermal fluctuations
%
% M = MOMENT(obj, 'option1', value1 ...)
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
%
% Options:
%
% nRand         The number of random Q points in the Brillouin-zone,
%               default is 1000.
% T             Temperature, default is zero.
% fid           Control the text output, default is one (Command Window).
% tol           Tolerance of the incommensurability of the magnetic
%               ordering wavevector. Deviations from integer values of the
%               ordering wavevector smaller than the tolerance are
%               considered to be commensurate. Default value is 1e-4.
% omega_tol     Tolerance on the energy difference of degenerate modes when
%               diagonalising the quadratic form, default is 1e-5.
%
% Output:
%
% M is a structure, with the following fields:
% moment        Size of the reduced moments, dimensions are [1 nMagExt].
% T           	Temperature.
% nRand         Number of random Q points.
% obj           The copy of the input obj.
%
% See also SW, SW.SPINWAVE, SW.GENMAGSTR.
%

inpForm.fname  = {'T'   'nRand' 'fid' 'tol' 'omega_tol'};
inpForm.defval = {0     1000    1     1e-4  1e-5       };
inpForm.size   = {[1 1] [1 1]   [1 1] [1 1] [1 1]      };

param = sw_readparam(inpForm, varargin{:});

% no modesorting
param.sortMode = false;

% size of the extended magnetic unit cell
nExt    = double(obj.mag_str.N_ext);
% magnetic ordering wavevector in the extended magnetic unit cell
km = obj.mag_str.k.*nExt;
% whether the structure is incommensurate
incomm = any(abs(km-round(km)) > param.tol);


fid      = param.fid;
nRand    = param.nRand;

% sum up the moment reduction
M.moment = 0;
hkl      = rand(3,nRand);

% Create the interaction matrix and atomic positions in the extended
% magnetic unit cell.
[SS, SI, ~] = obj.intmatrix;

% Introduce the opposite couplings.
% (i-->j) and (j-->i)
% transpose the JJ matrix as well [1 2 3 4 5 6 7 8 9] --> [6 9 12 7 10 13 8 11 14]
SS.new         = [SS.all(1:3,:)   -SS.all(1:3,:)  ];
SS.new(4:5,:)  = [SS.all([4 5],:)  SS.all([5 4],:)];
SS.new(6:14,:) = [SS.all(6:14,:)   SS.all([6 9 12 7 10 13 8 11 14],:) ]/2;
SS.all         = SS.new;

% Converts wavevctor list into the extended unit cell
hklExt  = bsxfun(@times,hkl,nExt')*2*pi;

% Calculates parameters eta and zed.
if isempty(obj.mag_str.S)
    error('sw:moment:NoMagneticStr','No magnetic structure defined in obj!');
end

M0 = obj.mag_str.S;
S0 = sqrt(sum(M0.^2,1));
% normal to rotation of the magnetic moments
n  = obj.mag_str.n;
nMagExt = size(M0,2);

if fid ~= 0
    if incomm
        fprintf(fid,'Calculating reduced moments of INCOMMENSURATE structure (nMagExt = %d, nRand = %d)...\n',nMagExt, nRand);
    else
        fprintf(fid,'Calculating reduced moments of COMMENSURATE structure (nMagExt = %d, nRand = %d)...\n',nMagExt, nRand);
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
MF  =  repmat(obj.unit.muB*SI.field*eta,[1 2]);

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

nSlice = ceil(nMagExt^2*nRand*6912/sw_freemem*2);

if nRand < nSlice
    if fid ~= 0
        fprintf(fid,'Memory allocation is not optimal, nMagExt is too large compared to the free memory!\n');
    end
    nSlice = nRand;
elseif nSlice > 1
    if fid ~= 0
        fprintf(fid,'To optimise memory allocation, Q is cut into %d pieces!\n',nSlice);
    end
end

hklIdx = [floor(((1:nSlice)-1)/nSlice*nRand)+1 nRand+1];

if fid == 1
    sw_status(0,1);
end

for jj = 1:nSlice
    % q indices selected for every chunk
    hklIdxMEM  = hklIdx(jj):(hklIdx(jj+1)-1);
    % q values contatining the k_m vector
    hklExtMEM  = hklExt(:,hklIdxMEM);
    nHklMEM = size(hklExtMEM,2);
    
    % Creates the matrix of exponential factors nCoupling x nHkl size.
    % Extends dR into 3 x 3 x nCoupling x nHkl
    ExpF = exp(1i*permute(sum(repmat(dR,[1 1 nHklMEM]).*repmat(permute(hklExtMEM,[1 3 2]),[1 nCoupling 1]),1),[2 3 1]))';
    
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
    
    g = diag([ones(nMagExt,1); -ones(nMagExt,1)]);
    
    % All the matrix calculations are according to White's paper
    
    gham = 0*ham;
    for ii = 1:nHklMEM
        gham(:,:,ii) = g*ham(:,:,ii);
    end
    
    [V, D] = eigorth(gham,param.omega_tol, param.sortMode);
    
    omega = zeros(2*nMagExt,nHklMEM);
    
    for ii = 1:nHklMEM
        % multiplication with g removed to get negative and positive
        % energies as well
        omega(:,ii) = D(:,ii);
        N           = diag(g*V(:,:,ii)'*g*V(:,:,ii));
        V(:,:,ii)   = V(:,:,ii)*diag(sqrt(1./N));
    end
    
    
    if param.T == 0
        nBose = double(omega<0);
    else
        nBose = double(omega>0).*(1./(exp(abs(omega)./(obj.unit.kB*param.T))-1) + 1);
    end
    
    % sums up the magnon population terms: U'*diag(n)*transpose(U)
    %M.moment = M.moment + sum(sum(conj(V).*repmat(permute(nBose,[3 1 2]),[nMagExt*2 1 1]).*V,3),2);
    M.moment = M.moment + sum(sum(abs(V).^2.*repmat(permute(nBose,[3 1 2]),[nMagExt*2 1 1]),3),2);
    
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% END MEMORY MANAGEMENT LOOP
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% sum up the result of <b^+b> and <bb^+> average over the first B.Z.
M.moment = S0 - ((M.moment(1:nMagExt)+M.moment(1+nMagExt:end))/nRand-1)'/2;
M.obj    = obj;
M.T      = param.T;
M.nRand  = nRand;

end