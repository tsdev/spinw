function chi = meanfield(obj, hkl, varargin)
% mean field calculation of the wave vector dependent susceptibility
%
% chi = MEANFIELD(obj, hkl, 'option1', value1, ...)
%
% Input:
%
% obj           Input structure, sw class object.
% hkl           Defines the Q points where chi is calculated, in reciprocal
%               lattice units, size is [3 nHkl]. Q can be also defined by
%               several linear scan in reciprocal space. In this case hkl
%               is cell type, where each element of the cell defines a
%               point in Q space. Linear scans are assumed between
%               consecutive points. Also the number of Q points can be
%               specified as a last element, it is 100 by defaults. For
%               example: hkl = {[0 0 0] [1 0 0]  50}, defines a scan along
%               (h,0,0) from 0 to 1 and 50 Q points are calculated along
%               the scan.
%
%               For symbolic calculation at a general reciprocal space
%               point use sym class input. For example to calculate chi
%               along (h,0,0): hkl = [sym('h') 0 0]. To do calculation at a
%               specific point do for example sym([0 1 0]), to calculate
%               the spectrum at (0,1,0).
%
% Options:
%
% Trel          Relative mean field temperature in Kelvin. An effective
%               temperature relative to the mean field critical temperature
%               Tc (the most negative eigenvalue of J(Q), the Fourier
%               transform of the exchange couplings). Default value is 0,
%               which means Tmf = Tc.
% Tc            Critical temperature, default value is calculated from the
%               exchange matrix sampled on the given Q points. If the Q
%               points don't contain the point where J(Q) is minimum, the
%               automatically determined Tc will be wrong. In this case it
%               is recommended to use this option to give the right Tc.
% formfact      Setting, that determines whether the magnetic form factor
%               is included in the spin-spin correlation function
%               calculation. Possible values:
%                   false   No magnetic form factor is applied (default).
%                   true    Magnetic form factors are applied, based on the
%                           label string of the magnetic ions, see sw_mff()
%                           function help.
%                   cell    Cell type that contains mixed labels and
%                           numbers for every symmetry inequivalent atom in
%                           the unit cell, the numbers are taken as
%                           constants.
%               For example: formfact = {0 'MCr3'}, this won't include
%               correlations on the first atom and using the form factor of
%               Cr3+ ion for the second atom.
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
%
% Output:
%
% chi           Structure with the following fields:
%
% Tc            Critical temperature in Kelvin, determined from the minimum
%               eigenvalue of J(Q) sampled on the Q points given in the
%               input.
% Tmf           Mean field temperature, Tmf = Trel + Tc.
% ...
%

inpForm.fname  = {'Trel' 'formfact' 'formfactfun' 'fitmode' 'gtensor' 'Tc' 'chi' };
inpForm.defval = {0      false       @sw_mff      false     false     []    []   };
inpForm.size   = {[1 1]  [1 -1]      [1 1]        [1 1]     [1 1]     [1 1] [1 1]};
inpForm.soft   = {false  false       false        false     false     true  true };

param = sw_readparam(inpForm, varargin{:});

fid = obj.fileid;

if ~isempty(param.chi)
    % use the previous calculation just change Tmf
    chi = param.chi;
    
    nHkl = size(chi.hkl,2);
    
    if fid ~= 0
        fprintf0(fid,['Realculating the mean field susceptibility for '...
            'different temperature using precalculated data (nHkl = %d)...\n'],nHkl);
    end

else
    
    % for linear scans create the Q line(s)
    if nargin > 1
        hkl = sw_qscan(hkl);
    else
        hkl = [];
    end
    
    % number of Q points
    nHkl = size(hkl,2);
    
    if fid ~= 0
        fprintf0(fid,['Calculating the mean field susceptibility '...
            '(nHkl = %d)...\n'],nHkl);
    end
    
    % calculate the Fourier transform of J
    chi = obj.fourier2(hkl);
    
    % calculates momentum transfer in A^-1 units.
    chi.hklA = 2*pi*(hkl'/obj.basisvector)';
    
    % solve the eigenvalue problem
    ffdim   = size(chi.ft);
    nMagExt = ffdim(3);
    ffdim1  = prod(ffdim([1 3]));
    % Hermitian matrix for the eigenvalue problem
    ff0 = reshape(permute(chi.ft,[1 3 2 4 5]),ffdim1,ffdim1,[]);
    
    % solve the iegenvalue problem
    [V, D] = eigorth(ff0);
    
    % energies (nMode x nHkl)
    chi.omega = real(D);
    
    % V: 3 x nMagExt x nMode x nHkl
    Sab = reshape(V,3,nMagExt,3*nMagExt,nHkl);
    
    % spin values are already included in the sw.fourier() function
    % % multiply the eigenmodes with the spin value
    % S   = obj.matom.S;
    % Sab = bsxfun(@times,V,S);
    
    % Create the interaction matrix and atomic positions in the extended
    % magnetic unit cell.
    if param.fitmode
        [~, SI] = obj.intmatrix('fitmode',2,'conjugate',true,'extend',false);
    else
        [~, SI] = obj.intmatrix('conjugate',true,'extend',false);
    end
    
    % message for magnetic form factor calculation
    if iscell(param.formfact) || param.formfact
        ffstrOut = 'The';
    else
        ffstrOut = 'No';
    end
    fprintf0(fid,[ffstrOut ' magnetic form factor is included in the spin-spin correlation function.\n']);
    
    % message for g-tensor calculation
    if param.gtensor
        gstrOut = 'The';
    else
        gstrOut = 'No';
    end
    fprintf0(fid,[gstrOut ' g-tensor is included in the spin-spin correlation function.\n']);
    
    % precalculate all magnetic form factors
    if iscell(param.formfact) || param.formfact
        if ~iscell(param.formfact)
            % unique atom labels
            uLabel = unique(obj.unit_cell.label(obj.unit_cell.S>0));
            % all atom labels
            aLabel = obj.unit_cell.label(obj.matom.idx);
            
            % save the form factor information in the output
            chi.formfact = obj.unit_cell.label(obj.unit_cell.S>0);
        else
            if numel(param.formfact) ~= unique(obj.matom.idx)
                error('sw:meanfield:WrongInput',['Number of form factor '...
                    'parameters has to equal to the number of symmetry inequivalent '...
                    'magnetic atoms in the unit cell!'])
            end
            % use the labels given as a cell input for all symmetry
            % inequivalent atom
            uLabel = param.formfact;
            aLabel = uLabel(obj.matom.idx);
            % convert numerical values to char() type
            aLabel = cellfun(@char,aLabel,'UniformOutput', false);
            
            % save the form factor information in the output
            chi.formfact = uLabel;
        end
        
        % stores the form factor values for each Q point and unique atom
        % label
        FF = zeros(nMagExt,nHkl);
        
        for ii = 1:numel(uLabel)
            lIdx = strcmp(aLabel,char(uLabel{ii}));
            FF(lIdx,:) = repmat(param.formfactfun(uLabel{ii},chi.hklA),[sum(lIdx) 1]);
        end
    else
        chi.formfact = false;
    end
    
    if param.gtensor
        % include the g-tensor
        gtensor = SI.g;
        Sab = permute(mmat(gtensor,permute(Sab,[1 5 2 3 4])),[1 3 4 5 2]);
    end
    
    if param.formfact
        % include the magnetic form factor
        Sab = bsxfun(@times,Sab,permute(FF,[3 1 4 2]));
    end
    chi.Sab = Sab;

end

if isempty(param.Tc)
    % find the critical temperature from the given Q points
    chi.Tc = sw_converter(2/3*abs(min(chi.omega(:))),'meV','K');
else
    chi.Tc = param.Tc;
end


% mean field temperature in Kelvin from the critical temperature and
% relative shift
chi.Tmf = param.Trel + chi.Tc; % K
Emf = sw_converter(chi.Tmf,'K','meV'); % meV


% denominator for the mean field susceptibility
deno = 3*Emf+2*chi.omega;

% wave vector dependent susceptibility
% 3 x 3 x nMagExt x nMagExt x nHkl
% assuming same moment sizes for all atoms!
chi.chi = sum(bsxfun(@rdivide,bsxfun(@times,permute(chi.Sab,[1 5 2 6 4 3]),...
    permute(conj(chi.Sab),[5 1 6 2 4 3])),permute(deno,[3 4 5 6 2 1])),6);

% neutron scattering cross section without plarisation factor
chi.int = real(squeeze(sumn(chi.chi,1:4)));

% divide Sab into symmetric and anti-symmetric components
chiS = (chi.chi + permute(chi.chi,[2 1 3 4 5]))/2;

% unplarised nuetron scattering cross section
% normalized scattering wavevector in xyz coordinate system.
hklAN = bsxfun(@rdivide,chi.hklA,sqrt(sum(chi.hklA.^2,1)));

% avoid NaN for Q=0
NaNidx = find(any(isnan(hklAN)));
for jj = 1:numel(NaNidx)
    if NaNidx(jj) < size(hklAN,2)
        hklAN(:,NaNidx(jj)) = hklAN(:,NaNidx(jj)+1);
    else
        hklAN(:,NaNidx(jj)) = [1;0;0];
    end
end

hkla = repmat(permute(hklAN,[1 3 4 5 2]),[1 3]);
hklb = repmat(permute(hklAN,[3 1 4 5 2]),[3 1]);

% Perpendicular part of the scattering wavevector.
qPerp = repmat(eye(3),[1 1 1 1 nHkl])- hkla.*hklb;

% Dynamical structure factor for neutron scattering
% Sperp: 1 x nHkl.
chi.Sperp = real(squeeze(sumn(bsxfun(@times,qPerp,chiS),1:4))');

if fid ~= 0
    fprintf0(fid,'Calculation finished.\n');
end


if ~param.fitmode
    chi.obj = copy(obj);
end

end