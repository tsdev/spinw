function spectra = swinc(obj, hkl, varargin)
% calculates spin-spin correlation function for planar incommensurate structures
%
% spectra = SWINC(obj, k, 'option1', value1 ...)
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
% n             Normal vector to the spin rotation plane, default is the
%               rotation plane of the spins, size is [1 3].
% sort          Whether to sort the modes using the function
%               eigenshuffle for plotting the dispersion and intensity.
%               Calculation time can be slightly reduced if sorting is
%               disabled. Default is true.
% formfact      Whether to include the magnetic form factor in the swConv
%               convoluted spectra. If true the form factor based on the
%               name of the first atom is used.
% fitmode       Speedup (for fitting mode only), default is false.
%
% Output:
%
% spectra is struct type, contains the calculated correlation function,
% with the following fields:
% omega     Calculated spin wave dispersion, size is [nMode nHkl],
%           where nMode==6*nMagExt for incommensurate structures, nMagExt
%           is the number of magnetic atoms in the extended unit cell. The
%           dispersion contains 2*nMagExt modes belonging to omega0(Q),
%           omega0(Q-km) and omega0(Q+km). Where omega0 is the calculated
%           dispersion. However the spin-spin correlation function for
%           neutron scattering contains the above three terms.
% Sab       Dynamical structure factor, dimensions are [3 3 nMode nHkl].
%           Each (:,:,i,j) submatrix contains the 9 correlation functions:
%           Sxx, Sxy, Sxz, ...
% hkl       Contains the input Q values, dimensions are [3 nHkl].
% hklA      Same Q values, but in reciproc Angstrom units in the lab
%           coordinate system, dimensions are [3 nHkl].
% ff        Magnetic form factor of the selected ion calculated for
%           hklA momentum transfer values, dimensions are [1 nHkl].
% n         Normal vector to the spin rotation plane, default is the
%           rotation plane of the spins, dimensions are [1 3].
% obj       The copy of the input obj.
%
% See also SW, SW.SPINWAVE, SW.POWSPEC, SW.OPTMAGSTR, SW_CONV.
%

% help when executed without argument
if nargin==1
    help sw.swinc
    return
end

% for linear scans create the Q line(s)
if iscell(hkl)
    hkl = sw_qscan(hkl);
end

nExt = double(obj.mag_str.N_ext);
% TODO
%kExt = mod(obj.mag_str.k.*nExt,1);
kExt = obj.mag_str.k.*nExt;

[n0, collinear] = sw_nvect(obj.mag_str.S);
if collinear
    n0 = obj.mag_str.n;
end

inpForm.fname  = {'n'   'sort' 'formfact' 'fitmode'};
inpForm.defval = {n0    true   true       false    };
inpForm.size   = {[1 3] [1 1]  [1 1]      [1 1]    };

param = sw_readparam(inpForm, varargin{:});

spectra   = struct;
nHkl      = size(hkl,2);

% Create the interaction matrix and atomic positions in the extended
% magnetic unit cell.
if param.fitmode
    [SS, SI, RR] = obj.intmatrix('fitmode',1);
else
    [SS, SI, RR] = obj.intmatrix;
end

% Converts wavevctor list into the extended unit cell
hklExt = bsxfun(@times,hkl,nExt')';

% Define transformation from the (n,u,v) coordinate system to (x,y,z)
% Descartes coordinate system.
tr      = zeros(3);
tr(:,1) = param.n';
[tr(:,2), tr(:,3)] = sw_cartesian(tr(:,1));

if isempty(obj.mag_str.S)
    error('sw:swinc:NoMagneticStr','No magnetic structure defined in swobj!');
end

% Calculates average spin value.
M0      = obj.mag_str.S;
nMagExt = size(M0,2);
avgS    = sum(sqrt(sum(M0.^2,1)))/nMagExt;


% Defines the spin angles in the (u,v) plane.
Mnuv = (M0'*tr)';
phi  = atan2(Mnuv(3,:),Mnuv(2,:));

% Transforms the single ion anisotropy matrices to the spin coordinate
% system.
SI.aninuv = SI.aniso*0;
for ii = 1:nMagExt
    [~, tr2] = sw_rot([1 0 0],phi(ii)-pi/2); % TODO
    SI.aninuv(:,:,ii) = tr2'* tr' * SI.aniso(:,:,ii) * tr * tr2;
end

% transforms the magnetic field into local coordinate system.
SI.effB = SI.field * tr * obj.unit.gamma;

% Construct a g matrix, x = (a,b,c,d,a+,b+,c+,d+) etc.
g = diag([ones(1,nMagExt) -ones(1,nMagExt)]);

% The matrix (size: 2*nMagExt x 2*nMagExt) stores the Hamiltonian for single
% k value, for the correlation functions the Hamiltonian is calculated for
% (k,k+km,k-km), there are nHkl number of different k vector along the scan
HAM  = zeros(nHkl,2*nMagExt,2*nMagExt,3);
gHAM = zeros(2*nMagExt,2*nMagExt,nHkl,3);

% Stores the eigenvectors and eigenvalues of HAM.
V     = zeros(2*nMagExt,2*nMagExt,nHkl,3);
% Dispersion contains all 2*nMagExt modes, for (k,k+km,k-km) vectors
omega = zeros(2*nMagExt,nHkl,3);

for Smode = 1:3
    % k values along the scan.
    k = hklExt + repmat((Smode-2)*kExt,nHkl,1);
    
    % Isortopic exchange interaction.
    % SSv contains a column of SS, Hiso(...) puts the 4x4 matrix into HAM
    for iConn = 1:size(SS.iso,2)
        SSv = SS.iso(:,iConn);
        HAM(:,:,:,Smode) = HisoV(SSv(6),HAM(:,:,:,Smode),kExt,nMagExt,k,RR(:,SSv(4))-RR(:,SSv(5))-SSv(1:3),phi(SSv(4))-phi(SSv(5)),SSv(4:5));
    end
    
    % Dzyaloshinskii-Moriya interaction.
    for iConn = 1:size(SS.dm,2)
        SSv = SS.dm(:,iConn);
        % D vector parallel with the normal of the spin-plane.
        D   = s.spin.n*SSv(6:8);
        HAM(:,:,:,Smode) = HdmV(D,HAM(:,:,:,Smode),kExt,nMagExt,k,RR(:,SSv(4))-RR(:,SSv(5))-SSv(1:3),phi(SSv(4))-phi(SSv(5)),SSv(4:5));
    end
    
    % Single-ion anisotropy energy.
    for ii = 1:nMagExt
        HAM(:,:,:,Smode) = Hsion(SI.aninuv(:,:,ii),HAM(:,:,:,Smode),nMagExt,ii);
    end
    
    % The S prefactor for the Hamiltonian, that is in all previous H functions.
    HAM(:,:,:,Smode) = avgS * HAM(:,:,:,Smode);
    
    % Zeeman energy in homogeneous magnetic field for commensurate structures.
    if any(SI.effB)
        warning('sw:swinc:FieldHelical','Magnetic field cannot be used for incommensurate helical structures!');
    end
    
    % Permutes back to the original order HAM dimensions.
    HAM = permute(HAM,[2 3 1 4]);
    
    % Matrix calculations according to White's paper.
    for n=1:nHkl
        gHAM(:,:,n,Smode) = g * HAM(:,:,n,Smode);
    end
    
    % Eigenshuffle sorts the eigenvalues and eigenvectors.
    if param.sort
        [Vtemp,Dtemp] = eigenshuffle(gHAM(:,:,:,Smode));
    else
        Vtemp = zeros(2*nMagExt,2*nMagExt,nHkl);
        Dtemp = zeros(2*nMagExt,nHkl);
        
        for ii = 1:nHkl
            [Vtemp(:,:,ii), D] = eig(gHAM(:,:,ii,Smode));
            Dtemp(:,ii)        = diag(D);
        end
    end
    
    V(:,:,:,Smode)     = Vtemp;
    omega(:,:,Smode)   = repmat(diag(g),1,nHkl).*Dtemp;
    
    for n = 1:nHkl
        M              = diag(g*V(:,:,n,Smode)'*g*V(:,:,n,Smode));
        V(:,:,n,Smode) = V(:,:,n,Smode)*diag(sqrt(1./M));
    end
    
    % Permutes HAM back.
    HAM = permute(HAM,[3 1 2 4]);
end

% Calculates the correlation functions.
Sn = [ones(1,nMagExt)  ones(1,nMagExt)];
Su = [    -sin(phi)    sin(phi)       ];
Sv = [     cos(phi)   -cos(phi)       ];

% Store the intensity of the transverse spin wave modes
Snn = zeros(2*nMagExt,nHkl);

% Incommensurate structures have Suv type correlations at wave vectors
% shifted by +/-k_m.
% TODO
% Su + Sv can give zero intensity !!!
Suv = zeros(2*nMagExt,nHkl,2);
for n=1:nHkl
    Suv(:,n,1) = (Su+Sv)/2*V(:,:,n,1).*conj((Su+Sv)/2*V(:,:,n,1));
    Suv(:,n,2) = (Su+Sv)/2*V(:,:,n,3).*conj((Su+Sv)/2*V(:,:,n,3));
    Snn(:,n)   = Sn*V(:,:,n,2).*conj(Sn*V(:,:,n,2));
end
% TODO 2x
Suv1 = avgS/2 * Suv(:,:,1)*2;
Suv2 = avgS/2 * Suv(:,:,2)*2;
Snn  = avgS/2 * Snn;

% use the rotating coordinate system
% SabRot = zeros(3,3,2*nMagExt,nHkl);
% SabRot(1,1,:,:) =  Snn;
% SabRot(2,2,:,:) =  (Suv1+Suv2)/2;
% SabRot(3,3,:,:) =  (Suv1+Suv2)/2;


% There are now three different Sab for k,k+/-km
% convert them: (nuv) --> (xyz) coordinates
% Snuv: 3 x nMode x nHkl x 3(q-km,q,q+km)
Snuv = zeros(3,2*nMagExt,nHkl,3);

Snuv(1,:,:,1) = Suv1;
Snuv(2,:,:,2) = Snn;
Snuv(3,:,:,3) = Suv2;

% Snuv: 3 x 3 x nMode x nHkl x 3(q-km,q,q+km)
Snuv = repmat(permute(Snuv,[1 5 2 3 4]),[1 3 1 1]);
% Sab0: 1 x 3 x nMode x nHkl x 3(q-km,q,q+km)
Sab0 = sum(repmat(inv(tr),[1 1 2*nMagExt nHkl 3]).*Snuv,1);
% Sab:  3 x 3 x nMode x nHkl x 3(q-km,q,q+km)
Sab = zeros(3,3,2*nMagExt,nHkl,3);
Sab(1,1,:,:,:) = Sab0(1,1,:,:,:);
Sab(2,2,:,:,:) = Sab0(1,2,:,:,:);
Sab(3,3,:,:,:) = Sab0(1,3,:,:,:);

% Scattering wave vector in lab coordinate system in Angstrom^-1 units
kA    = 2*pi*hkl'/obj.basisvector;

% Transforms kA vectors to (n,u,v) coordinate system.
%kAnuv = kA*tr;

% Intensity calculated for single unit cell, reduction factor follows.
iFact = 1/prod(nExt);


% Creates output structure with the calculated values.
%spectra.omega  = omega;
spectra.omega  = [omega(:,:,1); omega(:,:,2); omega(:,:,3)];
%spectra.SabRot = SabRot*iFact;
spectra.Sab    = cat(3,Sab(:,:,:,:,1),Sab(:,:,:,:,2),Sab(:,:,:,:,3))*iFact;
spectra.hkl    = hkl;
spectra.hklA   = kA';

if ~param.fitmode
    spectra.param  = param;
    spectra.obj    = copy(obj);
end

end


function [HAM] = HisoV(J,HAM,km,nMagExt,k,dr,dphi,Sconnect)
% function Hiso(km,nMagExt,k,dr,dphi,Sconnect) calculates the matrix form of the Heiseberg Hamiltonian of
% 2 spin scalar product
% km:   ordering wave vector
% nMagExt: number of sublattices
% k:    momentum transfer in Angstrom^-1
% dr:   the distance of the 2 spins
% phi:  the angle between the 2 spins - km*dr !!!
% Sconnect=[i j] contains the index of the two interacting spins
%

M1=J*cos( km*dr*2*pi+dphi);
M2=J*cos((km*dr*2*pi+dphi)/2)^2.* exp(1i*k*dr*2*pi);
M3=J*sin((km*dr*2*pi+dphi)/2)^2.* exp(1i*k*dr*2*pi);

i1=Sconnect(1);
i2=nMagExt+i1;
j1=Sconnect(2);
j2=nMagExt+j1;

HAM(:,i1,i1)=HAM(:,i1,i1)-M1;
HAM(:,i2,i2)=HAM(:,i2,i2)-M1;
HAM(:,j1,j1)=HAM(:,j1,j1)-M1;
HAM(:,j2,j2)=HAM(:,j2,j2)-M1;

HAM(:,i1,j1)=HAM(:,i1,j1)+M2;
HAM(:,i2,j2)=HAM(:,i2,j2)+M2;
HAM(:,j1,i1)=HAM(:,j1,i1)+conj(M2);
HAM(:,j2,i2)=HAM(:,j2,i2)+conj(M2);

HAM(:,i1,j2)=HAM(:,i1,j2)+M3;
HAM(:,i2,j1)=HAM(:,i2,j1)+M3;
HAM(:,j1,i2)=HAM(:,j1,i2)+conj(M3);
HAM(:,j2,i1)=HAM(:,j2,i1)+conj(M3);

end

% function [HAM] = Zeeman(Beff,HAM,nMagExt,phi,Sindex)
% % Zeeman function to calculate the Zeeman term
% %
%
% H1=cos(phi)*Beff(2)+sin(phi)*Beff(3);
%
% i1 = Sindex;
% i2 = nMagExt + i1;
%
% HAM(:,i1,i1) = HAM(:,i1,i1) + H1;
% HAM(:,i2,i2) = HAM(:,i2,i2) + H1;
%
% end


function [HAM] = Hsion(AA,HAM,nMagExt,Sindex)
%Hsion function to calculate the single ion anisotropy term
%

i1 = Sindex;
i2 = nMagExt + i1;

AA = AA/2;

A1 = AA(1,1) - AA(2,2) - 1i*(AA(1,2) + AA(2,1));
A2 = AA(1,1) - AA(2,2) + 1i*(AA(1,2) + AA(2,1));
A3 = AA(1,1) + AA(2,2) -  2* AA(3,3);

HAM(:,i1,i1) = HAM(:,i1,i1) + A3;
HAM(:,i2,i2) = HAM(:,i2,i2) + A3;

HAM(:,i1,i2) = HAM(:,i1,i2) + A1;
HAM(:,i2,i1) = HAM(:,i2,i1) + A2;

end

function [HAM]=HdmV(D,HAM,km,nMagExt,k,dr,dphi,Sconnect)
% function Hiso(km,nMagExt,k,dr,dphi,Sconnect) calculates the matrix form of
% the DM interacion where D is perpendicular to the spin plane of
% 2 spin vector product x-component
% km:   ordering wave vector
% nMagExt: number of sublattices
% k:    momentum transfer in Angstrom^-1
% dr:   the distance of the 2 spins
% phi:  the angle between the 2 spins - km*dr !!!
% Sconnect=[i j] contains the index of the two interacting spins
%

D1 = D/2*sin(km*dr*2*pi+dphi);
D2 = D/4*sin(km*dr*2*pi+dphi).* exp(1i*k*dr*2*pi);

i1 = Sconnect(1);
i2 = nMagExt+i1;
j1 = Sconnect(2);
j2 = nMagExt+j1;

HAM(:,i1,i1)=HAM(:,i1,i1) - D1;
HAM(:,i2,i2)=HAM(:,i2,i2) - D1;
HAM(:,j1,j1)=HAM(:,j1,j1) - D1;
HAM(:,j2,j2)=HAM(:,j2,j2) - D1;


HAM(:,i1,j1) = HAM(:,i1,j1) + D2;
HAM(:,i2,j2) = HAM(:,i2,j2) + conj(D2);
HAM(:,j1,i1) = HAM(:,j1,i1) + conj(D2);
HAM(:,j2,i2) = HAM(:,j2,i2) + D2;

HAM(:,i1,j2) = HAM(:,i1,j2) - D2;
HAM(:,i2,j1) = HAM(:,i2,j1) - conj(D2);
HAM(:,j1,i2) = HAM(:,j1,i2) - D2;
HAM(:,j2,i1) = HAM(:,j2,i1) - conj(D2);

end