function chi = mf(obj, hkl, varargin)
% mean field calculation of the wave vector dependent susceptibility
%
% chi = MF(obj, hkl, 'option1', value1, ...)
%
% Options:
%
% Tmf           Mean field temperature in Kelvin. An effective temperature,
%               has to be larger than the most negative eigenvalue of J(Q),
%               the Fourier transform of the exchange couplings. Default
%               value is 0.
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

inpForm.fname  = {'Tmf' 'formfact' 'formfactfun'};
inpForm.defval = {0     false       @sw_mff     };
inpForm.size   = {[1 1] [1 -1]      [1 1]       };

param = sw_readparam(inpForm, varargin{:});

% calculate the Fourier transform of J
chi = obj.fourier2(hkl);

% solve the eigenvalue problem
ffdim   = size(chi.ft);
nMagExt = ffdim(3);
ffdim1  = prod(ffdim([1 3]));
% Hermitian matrix for the eigenvalue problem
ff0 = reshape(permute(chi.ft,[1 3 2 4 5]),ffdim1,ffdim1,[]);

[V, D] = eigorth(ff0);

V = reshape(V,3,nMagExt,3*nMagExt,ffdim(5));

% calculate the wave vector dependent susceptibility

% mean field temperature in Kelvin
Tmf = param.Tmf; % K
Emf = sw_converter(Tmf,'K','meV'); % meV

% denominator for the MF susceptibility
deno = 3*Emf+2*D;
% wave vector dependent susceptibility
% 3 x 3 x nMagExt x nMagExt x nHkl
% assuming same moment sizes for all atoms!
chi.chi = sum(bsxfun(@rdivide,bsxfun(@times,permute(V,[1 5 2 6 4 3]),...
    permute(conj(V),[5 1 6 2 4 3])),permute(deno,[3 4 5 6 2 1])),6);

% neutron scattering cross section without plarisation factor
chi.int = squeeze(sumn(chi.chi,1:4));


% message for magnetic form factor calculation
if iscell(param.formfact) || param.formfact
    ffstrOut = 'The';
else
    ffstrOut = 'No';
end
fprintf0(fid,[ffstrOut ' magnetic form factor is included in the spin-spin correlation function.\n']);

% precalculate all magnetic form factors
if iscell(param.formfact) || param.formfact
    if ~iscell(param.formfact)
        % unique atom labels
        uLabel = unique(obj.unit_cell.label(obj.unit_cell.S>0));
        % all atom labels
        aLabel = obj.unit_cell.label(obj.matom(param.fitmode).idx);
        
        % save the form factor information in the output
        chi.formfact = obj.unit_cell.label(obj.unit_cell.S>0);
    else
        if numel(param.formfact) ~= unique(obj.matom.idx)
            error('sw:mf:WrongInput',['Number of form factor '...
                'parameters has to equal to the number of symmetry inequivalent '...
                'magnetic atoms in the unit cell!'])
        end
        % use the labels given as a cell input for all symmetry
        % inequivalent atom
        uLabel = param.formfact;
        aLabel = uLabel(obj.matom(param.fitmode).idx);
        % convert numerical values to char() type
        aLabel = cellfun(@char,aLabel,'UniformOutput', false);
        
        % save the form factor information in the output
        spectra.formfact = uLabel;
    end
    
    % stores the form factor values for each Q point and unique atom
    % label
    FF = zeros(nMagExt,nHkl);
    % Angstrom^-1 units for Q
    hklA0 = 2*pi*(hkl0'/obj.basisvector)';
    
    for ii = 1:numel(uLabel)
        lIdx = repmat(strcmp(aLabel,char(uLabel{ii})),[1 prod(nExt)]);
        FF(lIdx,:) = repmat(param.formfactfun(uLabel{ii},hklA0),[sum(lIdx) 1]);
    end
else
    spectra.formfact = false;
end



if iscell(param.formfact) || param.formfact
    % include the form factor in the z^alpha, z^beta matrices
    zeda = zeda.*repmat(permute(FF(:,hklIdxMEM),[3 4 5 1 2]),[3 3 2*nMagExt 2 1]);
    zedb = zedb.*repmat(permute(FF(:,hklIdxMEM),[3 4 1 5 2]),[3 3 2 2*nMagExt 1]);
end


% Calculates momentum transfer in A^-1 units.
chi.hklA = 2*pi*(hkl'/obj.basisvector)';

% divide Sab into symmetric and anti-symmetric components
chiS = (chi.chi + permute(chi.chi,[2 1 3 4 5]))/2;

% number of Q vectors
nHkl = size(hkl,2);

% unplarised nuetron scattering
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
chi.Sperp = squeeze(sumn(bsxfun(@times,qPerp,chiS),1:4))';

% include magnetic form factor









end