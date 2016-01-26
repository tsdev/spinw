function sFact = structfact(obj, hkl, varargin)
% calculates magnetic structure factor using FFT
%
% sFact = STRUCTFACT(obj, k, option1, value1, ...)
%
% Input:
%
% obj       Input sw object, contains positions of the magnetic atoms,
%           size of the magnetic supercell and the vector components of the
%           spins anf g-tensors.
% hkl       Defines the reciprocal lattice vectors where the magnetic
%           intensity is calculated. For commensurate structures these are
%           the possible positions of the magnetic Bragg peaks. For
%           incommensurate helical/conical structures 3 Bragg peaks
%           positions are possible: (k-km,k,k+km) around every reciprocal
%           lattice vector. In this case still the integer positions have
%           to be given and the code calculates the intensities at all
%           three points.
%
% Options:
%
% gtensor       If true, the g-tensor will be included in the static spin
%               correlation function. Including anisotropic g-tensor or
%               different g-tensor for different ions is only possible here.
%               Including a simple isotropic g-tensor is possible afterwards
%               using the sw_instrument() function.
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
% tol           Tolerance of the incommensurability of the magnetic
%               ordering wavevector. Deviations from integer values of the
%               ordering wavevector smaller than the tolerance are considered
%               to be commensurate. Default value is 1e-4.
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
% See also SW_PLOTSF, SW_INTSF, SW.ANNEAL, SW.GENMAGSTR.
%


inpF.fname  = {'gtensor' 'tol' 'formfact' 'formfactfun'};
inpF.defval = {false     1e-4  false       @sw_mff     };
inpF.size   = {[1 1]     [1 1] [1 -1]      [1 1]       };

param = sw_readparam(inpF, varargin{:});

sGrid = size(hkl);

matom = obj.matom;
nExt  = double(obj.mag_str.N_ext);
km    = obj.mag_str.k;
n     = obj.mag_str.n;

% convert km into the lu of the magnetic supercell
kmExt = mod(km.*nExt,1);

if sGrid(1) ~=3
    error('sw:structfact:WrongInput','The k input has to have 3 elements alogn the first dimension!')
end

if numel(sGrid)>2 && sGrid(1)==3
    % change grid into a list of k-points, data will be changed back in the
    % end.
    hkl = reshape(hkl,3,[]);
end

% check whether all hkl are integer in the supercell
hklExt = bsxfun(@times,hkl,nExt');
if sum(abs(hklExt(:)-round(hklExt(:))))/numel(hklExt) > 1e-10
    error('sw:structfact:WrongInput','All hkl values have to be integers.')
end

% whether the structure is incommensurate
incomm = any(abs(kmExt-round(kmExt)) > param.tol);

% special case if 2*km=tau structure
if incomm && all(abs(2*kmExt-round(2*kmExt))<param.tol)
    incomm = 2;
end

switch incomm
    case 0
        hklExtQ = hklExt;
        hklQ    = hkl;
    case 1
        % create (km-Q,km,km+Q) shifts for incommensurate structures
        hklExtQ = [bsxfun(@minus,hklExt,kmExt') hklExt bsxfun(@plus,hklExt,kmExt')];
        hklQ    = [bsxfun(@minus,hkl,(km./nExt)') hkl bsxfun(@plus,hkl,(kmExt./nExt)')];
    case 2
        % create (km+2*Q,km,km+Q) shifts for incommensurate structures
        hklExtQ = [bsxfun(@plus,hklExt,2*kmExt') hklExt bsxfun(@plus,hklExt,kmExt')];
        hklQ    = [bsxfun(@plus,hkl,2*(km./nExt)') hkl bsxfun(@plus,hklExt,(km./nExt)')];
end

% positions of the magnetic atoms in the supercell in lu units
matomExt = sw_extendlattice(nExt, matom);
% positions in unit of the crystallographic unit cell
RRext = matomExt.RRext;

% spins
S = obj.mag_str.S;

% calculate magnetic moments if gtensor is true, default g-value is 2
if param.gtensor
    % M = g*S
    [~,SI] = obj.intmatrix;
    %     mat = SI.;
    %     gIdx = double(obj.single_ion.g);
    %     gIdx(gIdx==0) = size(mat,3);
    M = permute(mmat(SI.g,permute(S,[1 3 2])),[1 3 2]);
else
    M = S;
end

% different dimensions
nHkl    = size(hklQ,2); 
nMagExt = obj.nmagext;

% Calculates momentum transfer in A^-1 units.
hklA = 2*pi*(hklQ'/obj.basisvector)';

% complex phases 1 x nHkl x nMagExt matrix
expF = exp(2*pi*1i*sum(bsxfun(@times,permute(hklExtQ,[3 2 4 1]),permute(RRext,[3 4 2 1])),4));

% calculate all magnetic form factors
if iscell(param.formfact) || param.formfact
    if ~iscell(param.formfact)
        % unique atom labels
        uLabel = unique(obj.unit_cell.label(obj.unit_cell.S>0));
        % all atom labels
        aLabel = obj.unit_cell.label(obj.matom.idx);
        
        % save the form factor information in the output
        sFact.formfact = obj.unit_cell.label(obj.unit_cell.S>0);
    else
        if numel(param.formfact) ~= numel(matom.idx)
            error('sw:spinwave:WrongInput',['Number of form factor '...
                'parameters has to equal to the number of magnetic '...
                'atoms in the unit cell!'])
        end
        % use the labels given as a cell input for all symmetry
        % inequivalent atom
        uLabel = param.formfact;
        aLabel = uLabel(obj.matom.idx);
        % convert numerical values to char() type
        aLabel = cellfun(@char,aLabel,'UniformOutput', false);
        
        % save the form factor information in the output
        sFact.formfact = uLabel;
    end
    
    % stores the form factor values for each Q point and unique atom
    % label
    FF = zeros(nMagExt,nHkl);
    
    for ii = 1:numel(uLabel)
        lIdx = repmat(strcmp(aLabel,char(uLabel{ii})),[1 prod(nExt)]);
        FF(lIdx,:) = repmat(param.formfactfun(uLabel{ii},hklA),[sum(lIdx) 1]);
    end
    
    % include magnetic form factor in the exponential prefactor
    expF = expF.*permute(FF,[3 2 1]);
    
else
    sFact.formfact = [];
end

% sum up spins
Mk = sum(bsxfun(@times,permute(M,[1 3 2]),expF),3);

if incomm>0
    % include rotation matrices
    nx  = [0 -n(3) n(2); n(3) 0 -n(1); -n(2) n(1) 0];
    nxn = n'*n;
    K1 = 1/2*(eye(3) - nxn - 1i*nx);
    K2 = nxn;
    
    Mk = [permute(sum(bsxfun(@times,K1,permute(Mk(:,1:end/3),[1 3 2])),2),[1 3 2]) ...
        permute(sum(bsxfun(@times,K2,permute(Mk(:,(end/3+1):end/3*2),[1 3 2])),2),[1 3 2]) ...
        permute(sum(bsxfun(@times,conj(K1),permute(Mk(:,(end/3*2+1):end),[1 3 2])),2),[1 3 2])];
end

if incomm == 2
    % sum up the km+Q and km+2Q values
    Mk = [Mk(:,(end/3+1):end/3*2) Mk(:,1:end/3)+Mk(:,(end/3*2+1):end)];
    
    hklQ = [hklQ(:,(end/3+1):end/3*2) hklQ(:,(end/3*2+1):end)];
end

% convert intensity per unit cell
Mk = Mk/prod(nExt);

% save output into struct
% structure factor
sFact.F2 = Mk.*conj(Mk);
% Fourier transform
sFact.Mk = Mk;
% momentum in rlu
sFact.hkl = hklQ;

sFact.hklA = 2*pi*(hklQ'/obj.basisvector)';

sFact.obj = copy(obj);
sFact.incomm = logical(incomm);

end