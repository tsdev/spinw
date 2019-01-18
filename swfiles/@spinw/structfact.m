function sFact = structfact(obj, kGrid, varargin)
% calculates magnetic and nuclear structure factor
% 
% ### Syntax
% 
% `sFact   = structfact(obj, kGrid,Name,Value)`
%
% `sfTable = structfact(obj, kGrid,Name,Value)`
%
% ### Description
% 
% `sFact   = structfact(obj, kGrid,Name,Value)` returns the calculated
% structure factors in units of barn. Magnetic structures (FM, AFM and
% helical) are checked against
% [FullProf](https://www.ill.eu/sites/fullprof/). The structure factor
% includes the site occupancy and Debye-Waller factors calculated from
% `obj.unit_cell.biso`, using the same definition as in FullProf.
% 
% ### Input Arguments
% 
% `obj`
% : [spinw] object.
% 
% `kGrid`
% : Defines the reciprocal lattice vectors where the structure
%      factor is to be calculated. For commensurate structures these
%      are the possible positions of the magnetic Bragg peaks. For
%      incommensurate helical/conical structures 3 Bragg peaks
%      positions are possible: $(\mathbf{k}-\mathbf{k}_m,\mathbf{k},\mathbf{k}+\mathbf{k}_m) around every reciprocal
%      lattice vector. In this case still the integer positions have
%      to be given and the code calculates the intensities at all
%      three points.
% 
% ### Name-Value Pair Arguments
% 
% `'mode'`
% : String, defines the type of calculation:
%   * `mag`     Magnetic structure factor and intensities for
%               unpolarised neutron scattering.
%   * `nucn`    Nuclear structure factor and neutron scattering
%               intensities.
%   * `nucx`    X-ray scattering structure factor and
%               intensities.
% 
% `'sortq'`
% : Sorting the reflections according to increasing momentum
%   value if `true`. Default is `false`.
% 
% `'formfact'`
% : If true, the magnetic form factor is included in the structure factor
%   calculation. The form factor coefficients are stored in
%   `obj.unit_cell.ff(1,:,atomIndex)`. Default value is `false`.
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
% `'gtensor'`
% : If true, the g-tensor will be included in the structure factor
%   calculation. Including anisotropic g-tensor or different
%   g-tensor for different ions is only possible here.
%
% `'lambda'`
% : Wavelength. If given, the $2\theta$ value for each reflection
%   is calculated.
% 
% `'dmin'`
% : Minimum $d$-value of a reflection, all higher order
%   reflections will be removed from the results.
% 
% `'output'`
% : String, defines the type of the output:
%   * `struct`  Results are returned in a struct type variable,
%               default.
%   * `table`   Results are returned in a table type output for
%               easy viewing and exporting.
% 
% `'tol'`
% : Tolerance of the incommensurability of the magnetic
%   ordering wavevector. Deviations from integer values of the
%   ordering wavevector smaller than the tolerance are considered
%   to be commensurate. Default value is $10^{-4}$.
% 
% `'fitmode'`
% : Speed up the calculation for fitting mode (omitting
%   cloning the [spinw] object into the output). Default is `false`.
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
% `sFact`
 % : Structure with the following fields:
%    * `F2`     Magnetic structure factor in a matrix with dimensions
%               $[3\times n_{hkl}]$.
%    * `Mk`     Square of the 3 dimensional magnetic structure factor,
%               dimensions are:
%               $[n_{ext}(1)\cdot f_{ext}(1)\times n_{ext}(2)\cdot f_{ext}(2)\times n_{ext}(3)\cdot f_{ext}(3)]$,
%               where $n_{ext}$ is the size of the extended unit cell.
%    * `hkl`    Contains the input $Q$ values in a matrix with dimensins of $[3\times n_{hkl}]$.
%    * `hklA`   Same as `hkl`, but in \\ang$^{-1}$ units in the
%               $xyz$ Cartesian coordinate system.
%    * `incomm` Whether the spectra calculated is incommensurate or not.
%    * `formfact` Cell containing the labels of the magnetic ions if form
%               factor in included in the spin-spin correlation function.
%    * `{tth}`  $2\theta$ value of the reflection for the given wavelength,
%               only given if a wavelength is provided.
%    * `obj`    Clone of the input `obj` object.
%
% `sfTable`
% : Table, optional output for quick viewing and saving the output into a
%   text file.
% 
% ### See Also
% 
% [sw_qgrid] \| [sw_plotsf] \| [sw_intsf] \| [spinw.anneal] \| [spinw.genmagstr]
%

inpF.fname  = {'mode' 'sortq' 'gtensor' 'formfact' 'formfactfun' 'fitmode'};
inpF.defval = {'mag'  false   false     false       @sw_mff      false    };
inpF.size   = {[1 -1] [1 1]   [1 1]     [1 1]       [1 1]        [1 1]    };
inpF.soft   = {false  false   false     false       false        false    };

inpF.fname  = [inpF.fname  {'lambda' 'output' 'dmin' 'rmzero' 'delta' 'fid'}];
inpF.defval = [inpF.defval {[]       'struct' []     false    1e-10   -1   }];
inpF.size   = [inpF.size   {[1 1]    [1 -2]   [1 1]  [1 1]    [1 1]   [1 1]}];
inpF.soft   = [inpF.soft   {true     false    true   false    false   false}];
  
param = sw_readparam(inpF, varargin{:});
pref = swpref;

if param.fid == -1
    fid = pref.fid;
else
    fid = param.fid;
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
        Wd = bsxfun(@times,biso(atom.idx),permute(1./(sFact.d.^2),[1 3 2]))/4;
        
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
    qA = sqrt(sum(sFact.hklA.^2,1));
    
    [qA,idx] = sort(qA,2);
    if ~isempty(param.dmin)
        idx = idx(2*pi./qA>param.dmin);
    end
    
    sFact.hklA = sFact.hklA(:,idx);
    sFact.d    = sFact.d(:,idx);
    
    if isfield(sFact,'Sab')
        sFact.Sab = reshape(sFact.Sab,3,3,[]);
        sFact.Sab = sFact.Sab(:,:,idx);
    end
    
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
    
    if param.rmzero
        idx = sFact.Sperp > param.delta;
        sFact.hkl  = sFact.hkl(:,idx);
        sFact.d    = sFact.d(:,idx);
        sFact.hklA = sFact.hklA(:,idx);
        sFact.Sperp = sFact.Sperp(1,idx);
        if isfield(sFact,'Sab')
            sFact.Sab = sFact.Sab(:,:,idx);
        end
        if isfield(sFact,'km')
            sFact.km = sFact.km(1,idx);
        end
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

if ~isempty(param.lambda)
    sFact.tth = 2*asind(sqrt(sum(sFact.hklA.^2,1))*param.lambda/4/pi);
end

% create output parameters
sFact.param    = param;
sFact.unit     = 'barn';
if ~param.fitmode
    sFact.obj = copy(obj);
end

% remove nan for Q=0
if isfield(sFact,'Sab')
    sFact.Sab(isnan(sFact.Sab)) = 0;
end
sFact.Sperp(isnan(sFact.Sperp)) = 0;

% create table if requested
switch param.output
    case 'table'
        sTab = table;
        sTab.h     = sFact.hkl(1,:)';
        sTab.k     = sFact.hkl(2,:)';
        sTab.l     = sFact.hkl(3,:)';
        switch param.mode
            case 'mag'
                dStr = 'magnetic neutron';
                sTab.km    = sFact.km';
            case 'nucn'
                dStr = 'nuclear neutron';
            case 'nucx'
                dStr = 'X-ray';
        end
        sTab.F2 = sFact.Sperp';
        sTab.d  = sFact.d';
        sTab.Properties.VariableUnits(1:3)  = {'r.l.u.' 'r.l.u.' 'r.l.u.'};
        sTab.Properties.VariableUnits{'F2'} = 'barn';
        sTab.Properties.VariableUnits{'d'}  = obj.unit.label{1};

        sTab.Properties.Description = ['Calculated ' dStr ' scattering cross section'];
        sTab.Properties.UserData.obj   = sFact.obj;
        sTab.Properties.UserData.param = sFact.param;
        if ~isempty(param.lambda)
            sTab.tth = sFact.tth';
            sTab.Properties.VariableUnits{'tth'} = 'degree';
        end
        
        sFact = sTab;
    case 'struct'
    otherwise
        error('spinw:structfact:WrongInput','''output'' option should be ''struct'' or ''table''!');
end

end