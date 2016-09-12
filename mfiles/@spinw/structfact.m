function sFact = structfact(obj, kGrid, varargin)
% calculates magnetic and nuclear structure factor
%
% sFact = STRUCTFACT(obj, kGrid, option1, value1, ...)
%
% The calculated structure factors are in barn units. Magnetic structures
% (FM, AFM and HELical) are checked against FullProf. The structure factor
% includes the site occupancy and Debye-Waller factors calculated from
% obj.unit_cell.biso, using the same definition as FullProf.
%
% Input:
%
% obj       Input spinw object, contains crystal and/or magnetic structure.
% kGrid     Defines the reciprocal lattice vectors where the structure
%           factor is to be calculated. For commensurate structures these
%           are the possible positions of the magnetic Bragg peaks. For
%           incommensurate helical/conical structures 3 Bragg peaks
%           positions are possible: (k-km,k,k+km) around every reciprocal
%           lattice vector. In this case still the integer positions have
%           to be given and the code calculates the intensities at all
%           three points.
%
% Options:
%
% mode          String, defines the type of calculation:
%                   mag     Magnetic structure factor and intensities for
%                           unpolarised neutron scattering.
%                   nucn    Nuclear structure factor and neutron scattering
%                           intensities.
%                   nucx    X-ray scattering structure factor and
%                           intensities.
% sortq         Sorting the reflections according to increasing momentum
%               value if true. Default is false.
%
% gtensor       If true, the g-tensor will be included in the static spin
%               correlation function, including anisotropic g-tensor or
%               different g-tensor per ion.
%
% formfact      If true, the magnetic form factor is included in the
%               spin-spin correlation function calculation. The form factor
%               coefficients are stored in obj.unit_cell.ff(1,:,atomIndex).
%               Default value is false.
%
% formfactfun   Function that calculates the magnetic form factor for given
%               Q value. Default value is @sw_mff(), that uses a tabulated
%               coefficients for the form factor calculation. For
%               anisotropic form factors a user defined function can be
%               written that has the following header:
%                   F = @formfactfun(atomLabel,Q)
%               where the parameters are:
%                   F   row vector containing the form factor for every
%                       input Q value
%                   atomLabel string, label of the selected magnetic atom
%                   Q   matrix with dimensions of [3 nQ], where each column
%                       contains a Q vector in Angstrom^-1 units.
% tol           Tolerance of the incommensurability of the magnetic
%               ordering wavevector. Deviations from integer values of the
%               ordering wavevector smaller than the tolerance are considered
%               to be commensurate. Default value is 1e-4.
% fitmode       Speed up the calculation for fitting mode (omitting
%               copying the spinw object to the output). Default is false.
%
% Output:
%
% 'sFact' is a structure with the following fields:
% F2            Magnetic structure factor in a matrix with dimensions
%               [3 x nHkl].
% Mk            Square of the 3 dimensional magnetic structure factor,
%               dimensions are:
%                  [nExt(1)*fExt(1) nExt(2)*fExt(2) nExt(3)*fExt(3)],
%               where nExt is the size of the extended unit cell.
% hkl           Contains the input Q values, dimensins are [3 nHkl].
% hklA          Same Q values, but in reciproc Angstrom units in the
%               lab coordinate system, dimensins are [3 nHkl].
% incomm        Whether the spectra calculated is incommensurate or not.
% formfact      Cell containing the labels of the magnetic ions if form
%               factor in included in the spin-spin correlation function.
% obj           Copy of the input obj object.
%
% See also SW_QGRID, SW_PLOTSF, SW_INTSF, SPINW.ANNEAL, SPINW.GENMAGSTR.
%

inpF.fname  = {'mode' 'sortq' 'gtensor' 'formfact' 'formfactfun' 'fitmode'};
inpF.defval = {'mag'  false   false     false       @sw_mff      false    };
inpF.size   = {[1 -1] [1 1]   [1 1]     [1 1]       [1 1]        [1 1]    };

param = sw_readparam(inpF, varargin{:});

if param.fitmode
    fid = 0;
else
    fid = obj.fileid;
end

% make a list of k-vectors from a grid and remember original dimensions
kDim  = size(kGrid);
hkl   = reshape(kGrid,3,[]);

% number of Q point
nQ    = size(hkl,2);

if kDim(1)~=3
    error('spinw:structfact:Wronginput','Dimensions of input hkl matrix are wrong!')
end

% constant for magnetic intensity
% neutron gyromagnetic ratio
gamma = 1.91304272;
% classical radius of the electron
r0 = 2.8179403267e-15; % m
% magnetic cross section constant in barn
constM = (gamma*r0/2)^2*1e28;

% occupancy
occ = obj.unit_cell.occ;
% isotropic displacement
biso = obj.unit_cell.biso;

switch param.mode
    case 'mag'
        % message for magnetic form factor calculation
        ffstrOut = {'No' 'The'};
        fprintf0(fid,[ffstrOut{param.formfact+1} ' magnetic form factor is'...
            ' included in the calculated structure factor.\n']);
        
        % magnetic atoms in the unit cell
        matom = obj.matom;
        % get k-vector
        km0  = obj.mag_str.k;
        % number of different magnetic propagation vectors
        nK0   = size(km0,2);
        % number of propagation vectors (count -km)
        nProp = any(mod(2*km0,1),1)+1;
        % number of prop vectors together with -km
        nK = sum(nProp);
        % store all km & -km
        kmStore = zeros(3,nK);
        iSign   = zeros(1,nK);
        iFact   = zeros(1,nK);
        kmIdx   = zeros(1,nK);
        
        % empty variables
        sFact.Sab   = zeros(3,3,nQ,nK);
        sFact.Sperp = zeros(nQ,nK);
        sFact.hklA  = zeros(3,nQ,nK);
        sFact.d     = zeros(1,nQ,nK);

        % loop over all prop vector to add the -km to the list
        idx = 1;
        for ii = 1:nK0
            if nProp(ii) == 1
                kmStore(:,idx) = km0(:,ii);
                iSign(1,idx)   = 1;
                iFact(1,idx)   = 2;
                kmIdx(1,idx)   = ii;
                idx = idx + 1;
            else
                idxL = idx+[0 1];
                kmStore(:,idxL) = [km0(:,ii) -km0(:,ii)];
                iSign(1,idxL)         = [1 -1];
                iFact(1,idxL)         = [1/2 1/2];
                kmIdx(1,idxL)         = [ii ii];
                idx = idx + 2;
            end
        end
        
        for ii = 1:nK
            % loop over all km vector one-by-one also including -km
            % select propagation vector
            km  = repmat(kmStore(:,ii),[1 nQ]);
            % position of magnetic reflections
            hklKm = bsxfun(@plus,hkl,km);
            % in Angstrom units [3 nQ nK]
            sFact.hklA(:,:,ii) = (hklKm'*obj.rl)';
            
            % Fourier components with the selected propagation vector
            % [3 nMagAtom 1]
            F1 = obj.mag_str.F(:,:,kmIdx(ii));
            % [3 nMagAtom] [3 _ nQ] --> F1*exp(i*kappa*r_ip) [3 nMagAtom nQ]
            kPerm = permute(iSign(ii)*hklKm,[1 3 2]);
            % include exp() factor
            F1 = bsxfun(@times,F1,exp(sum(bsxfun(@times,1i*2*pi*matom.r,kPerm),1)));

            % d-spacing in Angstrom units [1 nQ nK]
            sFact.d(:,:,ii) = 2*pi./sqrt(sum(sFact.hklA(:,:,ii).^2,1));
            % Debye-Waller factor [1 nMagAtom nQ]
            Wd = exp(-bsxfun(@times,biso(matom.idx),permute(1./(sFact.d(:,:,ii).^2),[1 3 2]))/4);
            
            % include occupancy and D-W factor
            F1 = bsxfun(@times,F1,bsxfun(@times,occ(matom.idx),Wd));
            
            if param.formfact
                % include magnetic form factors
                % store form factor per Q point for each atom in the
                % magnetic cell
                % [nMagAtom nQ]
                FF = param.formfactfun(permute(obj.unit_cell.ff(1,:,matom.idx),...
                    [3 2 1]),sFact.hklA(:,:,ii));
                % sum over magnetic atoms
                F1 = sum(bsxfun(@times,F1,permute(FF,[3 1 2])),2);
            else
                F1 = sum(F1,2);
            end
            
            % S^ab for elastic scattering [3 3 nQ]
            F2 = bsxfun(@times,permute(F1,[1 2 3]),conj(permute(F1,[2 1 3])));
            % [3 nQ]
            hklAnorm = bsxfun(@rdivide,sFact.hklA(:,:,ii),sqrt(sum(sFact.hklA(:,:,ii).^2,1)));
            % 1-q^2 [3 3 nQ]
            mq = bsxfun(@minus,eye(3),bsxfun(@times,permute(hklAnorm,[1 3 2]),permute(hklAnorm,[3 1 2])));
            % [3 3 nQ nK]
            sFact.Sab(:,:,:,ii) = iFact(ii)*F2.*mq;
            % sum up for non-polarized calculation, keep only the real part
            sFact.Sperp(:,ii) = constM*real(permute(sumn(sFact.Sab(:,:,:,ii),[1 2]),[3 1 2 4]));
        end

    case 'nucn'
        % nuclear structure factor
        % including occupancy & isotropic displacement

        % precalculation the atoms in the unit cell
        atom = obj.atom;
        % scattering length, convert to sqrt(barn)
        bc = obj.unit_cell.b(1,atom.idx)*0.1;
        
        sFact.hklA = (hkl'*obj.rl)';
        % d-spacing in Angstrom units [1 nQ]
        sFact.d = 2*pi./sqrt(sum(sFact.hklA.^2,1));
        % Debye-Waller factor [1 nAtom nQ]
        Wd = bsxfun(@times,biso(atom.idx),permute(1./(d.^2),[1 3 2]))/4;
        
        % nuclear unit-cell structure factor (fast, but takes lots of memory)
        % [1 1 nQ]
        F1 = sum(bsxfun(@times,bc.*occ(atom.idx),exp(-Wd+2*pi*1i*sum(bsxfun(@times,atom.r,permute(hkl,[1 3 2])),1))),2);
        % cross section
        sFact.Sperp = permute(F1.*conj(F1),[3 1 2]);
        
        % no magnetic wave vector
        nK = 1;
        % no magnetic form factor
        param.formfact = false;
        
    case 'nucx'
        fprintf('X-ray structure factor calculation is not implemented yet!')
    otherwise
        error('spinw:structfact:WrongInput','Wrong ''mode'' string!')
end

if param.sortq
    % don't resize the matrices, but sort
    sFact.hklA = reshape(sFact.hklA,3,[]);
    [~,idx] = sort(sum(sFact.hklA.^2,1));
    
    sFact.hklA = sFact.hklA(:,idx);
    
    sFact.Sab = reshape(sFact.Sab,3,3,[]);
    sFact.Sab = sFact.Sab(:,:,idx);
    
    sFact.Sperp = reshape(sFact.Sperp,1,[]);
    sFact.Sperp = sFact.Sperp(:,idx);
    % store reciprocal lattice vectors
    sFact.hkl = repmat(hkl,[1 nK]);
    sFact.hkl = sFact.hkl(:,idx);
    switch param.mode
        case 'mag'
            % store the magnetic propagation vector
            sFact.km  = repmat(permute(1:nK,[1 3 2]),[1 nQ 1]);
            sFact.km  = reshape(sFact.km,1,[]);
            sFact.km  = sFact.km(:,idx);
    end
else
    % resize the matrices to the dimensions of the input q-grid
    % [nQ1 nQ2 nQ3 nK]
    sFact.Sperp = reshape(sFact.Sperp,[kDim(2:end) nK]);
    % store reciprocal lattice vectors
    sFact.hkl = kGrid;
    switch param.mode
        case 'mag'
            sFact.km    = kmStore;
            % [3 3 nQ1 nQ2 nQ3 nK]
            sFact.Sab = reshape(sFact.Sab,[3 3 kDim(2:end) nK]);
    end
    
end

% create output parameters
sFact.param    = param;
sFact.unit     = 'barn';
if ~param.fitmode
    sFact.obj      = copy(obj);
end

end