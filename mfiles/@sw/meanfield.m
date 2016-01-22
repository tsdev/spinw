function chi = meanfield(obj, hkl, varargin)
% mean field calculation of the wave vector dependent susceptibility
%
% chi = MEANFIELD(obj, hkl, 'option1', value1, ...)
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

inpForm.fname  = {'Tmf' 'formfact' 'formfactfun' 'fitmode' 'gtensor'};
inpForm.defval = {0     false       @sw_mff      false     false    };
inpForm.size   = {[1 1] [1 -1]      [1 1]        [1 1]     [1 1]    };

param = sw_readparam(inpForm, varargin{:});

fid = obj.fileid;

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
hklA = 2*pi*(hkl'/obj.basisvector)';

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
        aLabel = obj.unit_cell.label(obj.matom(param.fitmode).idx);
        
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
        aLabel = uLabel(obj.matom(param.fitmode).idx);
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
        FF(lIdx,:) = repmat(param.formfactfun(uLabel{ii},hklA),[sum(lIdx) 1]);
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

% mean field temperature in Kelvin
Tmf = param.Tmf; % K
Emf = sw_converter(Tmf,'K','meV'); % meV

% denominator for the mean field susceptibility
% correction, no factor-2 in front of D
deno = 3*Emf+D;

% wave vector dependent susceptibility
% 3 x 3 x nMagExt x nMagExt x nHkl
% assuming same moment sizes for all atoms!
chi.chi = sum(bsxfun(@rdivide,bsxfun(@times,permute(Sab,[1 5 2 6 4 3]),...
    permute(conj(Sab),[5 1 6 2 4 3])),permute(deno,[3 4 5 6 2 1])),6);

% neutron scattering cross section without plarisation factor
chi.int = real(squeeze(sumn(chi.chi,1:4)));

% divide Sab into symmetric and anti-symmetric components
chiS = (chi.chi + permute(chi.chi,[2 1 3 4 5]))/2;

% unplarised nuetron scattering cross section
% normalized scattering wavevector in xyz coordinate system.
hklAN = bsxfun(@rdivide,hklA,sqrt(sum(hklA.^2,1)));

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

% save variables
chi.hklA = hklA;

if ~param.fitmode
    chi.obj = copy(obj);
end

end