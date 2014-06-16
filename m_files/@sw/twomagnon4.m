function spectra = twomagnon(obj, hkl, varargin)
% calculates two magnon scattering using linear spin wave theory
%
% spectra = TWOMAGNON(obj, k, 'option1', value1 ...)
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
% nRand         Number of Monte Carlo cycles per Q points. Default is 1e4.
% T             Temperature, for calculating the Bose distribution of
%               magnons for the spin-spin correlation function. Default is
%               taken from obj.single_ion.T value.
% Evect         Energy grid for binning the correlation function along
%               energy. Default is linspace(0,1,101);
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
%
% Output:
%
% spectra is a structure, with the following fields:
% omega         Calculated spin wave dispersion, dimensins are
%               [nMode nHkl], where nMagExt is the number of magnetic
%               atoms in the extended unit cell.
% V             Transformation matrix between the magnon operators and the
%               normal magnon operators. Dimensions are [nMode nMode nHkl].
%
% nMode is the number of magnetic mode, it is double the number of magnetic
% atoms in the magnetic cell/supercell.
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
    help sw.twomagnon
    return
end


% for linear scans create the Q line(s)
if iscell(hkl)
    hkl = sw_qscan(hkl);
elseif numel(hkl)==3
    hkl = hkl(:);
end

T0 = obj.single_ion.T;
E0 = linspace(0,1,101);

inpForm.fname  = {'nRand' 'T'  'Evect' 'fid' 'tol'  'omega_tol' 'hermit' };
inpForm.defval = {1e4     T0    E0      1      1e-4  1e-5        true    };
inpForm.size   = {[1 1]   [1 1] [1 -1]  [1 1]  [1 1] [1 1]       [1 1]   };

param = sw_readparam(inpForm, varargin{:});

fid = param.fid;

% size of the extended magnetic unit cell
nExt    = double(obj.mag_str.N_ext);
% magnetic ordering wavevector in the extended magnetic unit cell
km = obj.mag_str.k.*nExt;
% normal to the spin rotation plane
n  = obj.mag_str.n;

% whether the structure is incommensurate
incomm = any(abs(km-round(km)) > param.tol);

if incomm
    tol = param.tol*2;
    if sum(abs(mod(abs(2*km)+tol,1)-tol).^2) < tol
        error('sw:spinwave:kHalf',['For magnetic structures where 2*km equal' ...
            ' to lattice translation vector, the magnetic unit cell has to'...
            ' be extended!']);
    end
    nKm = 3;
    kmList = [-km' zeros(3,1) km'];
else
    nKm = 1;
    kmList = zeros(3,1);
end

% number of Q points
nHkl = size(hkl,2);
nRand = param.nRand;

% bins along energy
Evect   = sort(param.Evect);
epsilon = 1e-8;

if ~isempty(Evect)
    Evect = [Evect(1)-epsilon; Evect(:); Evect(end)+epsilon];
else
    Evect = [-epsilon; epsilon];
end
nE      = numel(Evect);

% eta vectors
M0 = obj.mag_str.S;
eta = bsxfun(@rdivide,M0,sqrt(sum(M0.^2,1)));

% number of omega modes
nMode = size(M0,2)*2;

% temperature in energy units
kBT = obj.unit.kB*param.T;

Sab2 = zeros(nE,nHkl,3,3);

if incomm
    % Rodrigues' rotation formula.
    nx  = [0 -n(3) n(2); n(3) 0 -n(1); -n(2) n(1) 0];
    nxn = n'*n;
    K1 = 1/2*(eye(3) - nxn - 1i*nx);
    K2 = nxn;
    
    % keep the rotation invariant part of Sab
    nx  = [0 -n(3) n(2);n(3) 0 -n(1);-n(2) n(1) 0];
    nxn = n'*n;
    m1  = eye(3);
end

if fid == 1
    sw_status(0,1);
end

% size of the extended magnetic unit cell
nExt    = double(obj.mag_str.N_ext);

mAtom = sw_extendlattice(nExt, obj.matom);
RR = mAtom.RRext;

% Converts wavevctor list into the extended unit cell
hklExt  = bsxfun(@times,hkl,nExt')*2*pi;

% exp(-i*k*(Ri-Rj) factor
ExpF = exp(1i*sum(repmat(permute(hklExt,[1 3 2]),[1 nMode/2 1]).*repmat(RR,[1 1 nHkl]),1));
%ExpF = repmat(bsxfun(@times,permute(ExpF,[2 1 3]),conj(ExpF)),[2 2 1]);

for ii = 1:nHkl
    for jj = 1:nKm
        % random Q points from the first Brillouin-zone
        randQ = rand(3,nRand);
        
        % calculate dispersion, V matrix
        [omega1, V1] = obj.corespec(randQ);
        [omega2, V2] = obj.corespec(bsxfun(@plus,-randQ,hkl(:,ii)+kmList(:,jj)));
        
        omega1 = omega1(1:nMode/2,:);
        omega2 = omega2(1:nMode/2,:);
        
        % energy of the two magnon process = sum(omega)
        omega = omega1+omega2;
        
        % intensity
        % dimensions [nMode nMode nRand alpha beta]
        ExpF1 = bsxfun(@times,ExpF(:,:,ii),permute(eta,[3 2 4 1]));
        
        Int = mmat(bsxfun(@times,ExpF1,permute(conj(V1(1:nMode/2,1:nMode/2,:)),[2 1 3])),permute(V2(1:nMode/2,1:nMode/2,:),[2 1 3]));
        Int0 = bsxfun(@times,Int,permute(conj(Int),[1 2 3 5 4]));
        
        Int = zeros(nMode/2,nRand,3,3);
        for kk = 1:nMode/2
            Int(kk,:,:,:) = Int0(kk,kk,:,:,:);
        end
        
        % Calculate Bose temperature factor for magnons
        % dimensions [nMode nMode nRand]
        if param.T==0
            nBose = double(omega1>0).*double(omega2>0);
        else
            nBose = 1./((exp(abs(omega1)./kBT)-1)+double(omega1>0).*(exp(abs(omega2)./kBT)-1)+double(omega2>0));
        end
        
        % intensity = intensity * Bose
        % dimensions [alpha beta nMode nRand]
        Int = permute(bsxfun(@times,Int,nBose),[3 4 1 2]);
        
        [~, idxE] = min(abs(repmat(real(omega),[1 1 nE])-repmat(permute(Evect,[2 3 1]),[nMode/2 nRand 1])),[],3);
        
        if incomm
            % introduce the rotation matrices from the Rodriguez formula
            switch jj
                case 1
                    Int = mmat(Int,K1);
                case 2
                    Int = mmat(Int,K2);
                case 3
                    Int = mmat(Int,conj(K1));
            end
        end
        
        for kk = 1:9
            [ind1,ind2]=ind2sub([3 3],kk);
            DSF = squeeze(Int(ind1,ind2,:,:));
            Sab2(:,ii,ind1,ind2) = Sab2(:,ii,ind1,ind2) + accumarray(idxE(:),DSF(:),[nE 1]);
        end
        
        if fid == 1
            sw_status(ii/nHkl*100);
        end
    end
end

Sab2 = permute(Sab2(2:end-1,:,:,:),[3 4 1 2]);

if incomm
    % integrate out the reference coordinate system
    Sab2 = 1/2*Sab2 - 1/2*mmat(mmat(nx,Sab2),nx) + ...
        1/2*mmat(mmat(nxn-m1,Sab2),nxn) + 1/2*mmat(mmat(nxn,Sab2),2*nxn-m1);
end

if fid == 1
    sw_status(100,2);
else
    if fid ~= 0
        fprintf(fid,'Calculation finished.\n');
    end
end

% Calculates momentum transfer in A^-1 units.
hklA = 2*pi*(hkl'/obj.basisvector)';

% Creates output structure with the calculated values.
spectra.Sab2   = Sab2/nRand;
spectra.hkl    = hkl;
spectra.hklA   = hklA;
spectra.incomm = incomm;
spectra.Evect  = Evect(2:end-1);
spectra.norm   = false;

% clone sw object
spectra.obj = copy(obj);

% save the important parameters
spectra.param = param;

end