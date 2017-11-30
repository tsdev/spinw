function spectra = sw_xray(spectra, varargin)
% calculates x-ray scattering cross section
% 
% ### Syntax
% 
% `spectra = sw_xray(spectra,Name,Value)`
% 
% ### Description
% 
% `spectra = sw_xray(spectra,Name,Value)` calculates x-ray scattering
% intensity for non-resonant inelastic x-ray scattering on phonons.
%  
% 
% ### Input Arguments
% 
% `spectra`
% : Input structure that contains the displacement-displacement
%   correlation function.
% 
% ### Name-Value Pair Arguments
%
% `'formfact'`
% : If true, the magnetic form factor is included in the spin-spin
%   correlation function calculation. The form factor coefficients are
%   stored in `obj.unit_cell.ff(1,:,atomIndex)`. Default value is `false`.
%
% `'formfactfun'`
% : Function that calculates the magnetic form factor for given $Q$ value.
%   value. Default value is `@sw_mff`, that uses a tabulated coefficients
%   for the form factor calculation. For anisotropic form factors a user
%   defined function can be written that has the following header:
%   ```
%   F = formfactfun(atomLabel,Q)
%   ```
%   where the parameters are:
%   * `F`           row vector containing the form factor for every input 
%                   $Q$ value
%   * `atomLabel`   string, label of the selected magnetic atom
%   * `Q`           matrix with dimensions of $[3\times n_Q]$, where each
%                   column contains a $Q$ vector in $\\ang^{-1}$ units.
%
% `'fid'`
% : Defines whether to provide text output. The default value is determined
%   by the `fid` preference stored in [swpref]. The possible values are:
%   * `0`   No text output is generated.
%   * `1`   Text output in the MATLAB Command Window.
%   * `fid` File ID provided by the `fopen` command, the output is written
%           into the opened file stream.
%
% ### Output Arguments
% 
% `spectra`
% : Structure that is same as the input with the following additional
%   fields:
%   * `param`   Input parameters.
%   * `Sperp`   $S_\perp(i_{mode},\mathbf{Q})$ x-ray scattering cross
%               section, stored in a matrix with dimensions of
%               $[n_{mode n_{hkl}]$.
% 
% ### See Also
% 
% [spinw] \| [spinw.spinwave] \| [sw_neutron]
%

if nargin == 0
    swhelp sw_xray
    return
end

inpForm.fname  = {'formfact' 'formfactfun' 'fid'};
inpForm.defval = {true        @sw_cff      -1   };
inpForm.size   = {[1 -1]      [1 1]        [1 1]};

param = sw_readparam(inpForm,varargin{:});

% loop over the twins
if ~iscell(spectra.Sab)
    Sab = {spectra.Sab};
else
    Sab = spectra.Sab;
end

nTwin = numel(Sab);
spectra.Sperp = cell(1,nTwin);

% copy out SpinW object
obj = spectra.obj;
% number of Q points
%nHkl = size(spectra.omega,2);

% convert Q into A-1 units
% Angstrom^-1 units for Q
QA = 2*pi*(spectra.hkl'/obj.basisvector)';
% calculate intensity
% atoms in the structure
atom = obj.atom;
% number of atoms in the unit cell
%nAtom = numel(atom.idx);

% get X-ray form factors
% precalculate all magnetic form factors
% calculate all magnetic form factors
if param.formfact
    spectra.formfact = true;
    % store form factor per Q point for each atom in the magnetic supercell
    FF = param.formfactfun(permute(obj.unit_cell.ff(2,:,atom.idx),[3 2 1]),QA);
else
    spectra.formfact = false;
end



% if param.formfact
%     if ~iscell(param.formfact)
%         % unique atom labels
%         uLabel = obj.unit_cell.label;
%         % all atom labels
%         aLabel = obj.unit_cell.label(atom.idx);
%         
%         % save the form factor information in the output
%         spectra.formfact = obj.unit_cell.label;
%     else
%         if numel(param.formfact) ~= numel(obj.unit_cell.label)
%             error('x_ray:WrongInput',['Number of form factor '...
%                 'parameters has to equal to the number of symmetry inequivalent '...
%                 'magnetic atoms in the unit cell!'])
%         end
%         % use the labels given as a cell input for all symmetry
%         % inequivalent atom
%         uLabel = param.formfact;
%         aLabel = uLabel(obj.atom.idx);
%         % convert numerical values to char() type
%         aLabel = cellfun(@char,aLabel,'UniformOutput', false);
%         
%         % save the form factor information in the output
%         spectra.formfact = uLabel;
%     end
%     
%     % stores the form factor values for each Q point and unique atom
%     % label
%     FF = zeros(nAtom,nHkl);
%     
%     for ii = 1:numel(uLabel)
%         lIdx = strcmp(aLabel,char(uLabel{ii}));
%         FF(lIdx,:) = repmat(param.formfactfun(uLabel{ii},QA),[sum(lIdx) 1]);
%     end
% else
%     spectra.formfact = false;
% end

% get the mass of the atoms
aLabel = obj.unit_cell.label(atom.idx);
% get mass
M = zeros(1,numel(aLabel));
for ii = 1:numel(aLabel)
    M(ii) = sw_atomdata(aLabel{ii},'mass');
end

% apply the X-ray scattering formula
for ii = 1:nTwin
    % apply form factor
    temp = bsxfun(@times,Sab{ii},permute(FF,[3 1 4 2]));
    % apply Debye-Waller factor if given
    if isfield(spectra,'DW')
        fprintf0(spectra.obj.fileid,'Including the Debye-Waller factor in the X-ray scattering cross section!\n');
        % QT*W*Q values
        QWQ = mmat(mmat(permute(QA,[3 1 4 2]),spectra.DW),permute(QA,[1 3 4 2]));
        temp = bsxfun(@times,temp,exp(-permute(QWQ,[1 3 2 4])));
    else
        fprintf0(spectra.obj.fileid,'No Debye-Waller factor is included in the X-ray scattering cross section!\n');
    end
    % apply 1/sqrt(M) factor on u(q)
    temp = bsxfun(@rdivide,temp,sqrt(M));
    % include phase factor
    phi = sum(bsxfun(@times,obj.atom.r,permute(spectra.hkl,[1 3 4 2])),1);
    temp = bsxfun(@times,temp,exp(1i*2*pi*phi));
    % apply e*Q factor
    spectra.Sperp{ii} = squeeze((abs(sum(bsxfun(@times,sum(temp,2),permute(QA,[1 3 4 2])),1))).^2);
    % apply 1/omega factor
    spectra.Sperp{ii} = spectra.Sperp{ii}./spectra.omega;
end

% remove cell if there are no twins
if nTwin == 1
    spectra.Sperp = spectra.Sperp{1};
end

end