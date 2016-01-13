function [w, s] = horace(obj, qh, qk, ql, p_in)
% calculates spin wave dispersion/correlation functions to be called from Horace
%
% [w, s] = HORACE(obj, qh, qk, ql, p)
%
% The function produces spin wave dispersion and intensity for Horace
% (<a href=http://horace.isis.rl.ac.uk>http://horace.isis.rl.ac.uk</a>).
%
% Input:
%
% obj           Input sw object.
% qh, qk, ql    Reciprocal lattice components in reciprocal lattice units.
% p             Parameters, in the order defined by horace_setpar
%
% Output:
%
% w         Cell that contains the spin wave energies. Every cell elements
%           contains a vector of spin wave energies for the corresponding
%           input Q vector.
% s         Cell that contains the calculated element of the spin-spin
%           correlation function. Every cell element contains a vector of
%           intensities in the same order as the spin wave energies in w.
%
% Example:
%
% ...
% horace_on;
% d3dobj = d3d(cryst.abc,[0 1 0 0],[0,0.01,1],[0 0 1 0],[0,0.01,1],[0 0 0 1],[0,0.1,10]);
% d3dobj = disp2sqw_eval(d3dobj,@cryst.horace,{'component','Sperp'},0.1);
% plot(d3dobj);
%
% This example creates a d3d object, a square in (h,k,0) plane and in
% energy between 0 and 10 meV. Then calculates the inelastice neutron
% scattering intensity of the spin wave model stored in cryst and plots it
% using sliceomatic.
%
% See also SW, SW.SPINWAVE, SW.MATPARSER, SW.HORACE_SETPAR, SW_READPARAM.
%

if nargin <= 1
    help sw.horace;
    return;
end

%{
if ~isempty(param.param)
    % change matrix values for Horace data fitting
    if nargin(param.func) < 0
        param.func(varargin{:});
    elseif nargin(param.func) == 2
        param.func(param.param);
    else
        error('sw:horace:WrongInput','User defined function with incompatible header!');
    end
end

if param.norm && param.dE == 0
    error('sw:horace:WrongInput',['To convert spin wave intensity to mbarn/meV/cell/sr'...
        ' units, give the energy bin step.'])
end

% calculate spin wave spectrum
if nargin > 5
    % include the fitmode option to speed up calculation
    if numel(varargin) == 1
        varargin{1}.fitmode = 2;
        spectra = obj.spinwave([qh(:) qk(:) ql(:)]',varargin{1});
    else
        spectra = obj.spinwave([qh(:) qk(:) ql(:)]',varargin{:},'fitmode',true);
    end
else
    spectra = obj.spinwave([qh(:) qk(:) ql(:)]','fitmode',true);
end
warning(warnState);
%}

% Set default parameters if the user has not already set them...
if ~isfield(obj.matrix,'horace')
    obj.horace_setpar();
end
% Changes the parameters required.
obj.matrix.horace.func('param',p_in,'mat',obj.matrix.horace.mapping, ...
      'selector',obj.matrix.horace.selector,'init',obj.matrix.horace.init);
% Sets the other parameters from the stored field.
fname = {'component' 'norm' 'dE'  'tol' 'hermit' 'notwin'}; % 'optmem'
fname = [fname  {'formfact' 'formfactfun' 'gtensor' 'omega_tol' 'useMex'}];
pars = {};
for ifl = 1:numel(fname)
    pars = [pars {fname{ifl} obj.matrix.horace.(fname{ifl})}];
end
param = obj.matrix.horace;

% Determines the number of k-points
hkl = [qh(:) qk(:) ql(:)]';
nHkl = size(hkl,2);

% Determines the optimum number of chunks
nMagExt = size(obj.mag_str.S,2);
if obj.matrix.horace.optmem == 0
    freeMem = sw_freemem;
    if freeMem > 0
        nSlice = ceil(nMagExt^2*nHkl*4000/freeMem);
    end
else
    nSlice = obj.matrix.horace.optmem;
end
hklIdx = [floor(((1:nSlice)-1)/nSlice*nHkl)+1 nHkl+1];
nTwin = obj.matrix.horace.notwin;
if nTwin==0
    nTwin = 1;
end
w = cell(1,2*nMagExt*nTwin);
s = cell(1,2*nMagExt*nTwin);
for jj = 1:2*nMagExt
    w{jj} = zeros(nHkl,1);
    s{jj} = zeros(nHkl,1);
end
% Chunk the k-point list to save memory.
for iSlice = 1:nSlice
    hklIdxMEM  = hklIdx(iSlice):(hklIdx(iSlice+1)-1);
    nHklMEM = numel(hklIdxMEM);
    % calculate Sab
    spectra = obj.spinwave(hkl(:,hklIdxMEM),'fitmode',3,pars{:},'optmem',1);
    % calculate Sperp
    spectra = sw_neutron(spectra);
    nMode = size(spectra.omega,1);

    % parse the component string
    if iscell(param.component)
        nConv = numel(param.component);
        parsed = cell(1,nConv);
        for ii = 1:numel(param.component)
            parsed{ii} = sw_parstr(param.component{ii});
        end
    elseif isstruct(param.component)
        nConv  = 1;
        parsed = {param.component};
        param.component = {parsed{1}.string};
    else
        nConv = 1;
        parsed = {sw_parstr(param.component)};
        param.component = {param.component};
    end

    % pack all cross section into a cell for easier looping
    if iscell(spectra.omega)
        nTwin = numel(spectra.omega);
        omega = spectra.omega;
        Sab = spectra.Sab;
        Sperp = spectra.Sperp;

    else
        nTwin = 1;
        omega = {spectra.omega};
        Sab = {spectra.Sab};
        Sperp = {spectra.Sperp};
    end

    % extract the requested cross section
    % DSF stores the intensity that is going to be convoluted
    DSF = cell(nConv,nTwin);
    % select value from Sperp or Sab matrices
    for tt = 1:nTwin
        for ii = 1:numel(parsed)
            par0 = parsed{ii};
            DSF{ii,tt} = zeros(nMode, nHklMEM);
            for jj = 1:length(par0.type)
                switch par0.type{jj}(1)
                    case 1
                        DSF{ii,tt} = DSF{ii,tt} + par0.preFact(jj)*Sperp{tt};
                    case 2
                        DSF{ii,tt} = DSF{ii,tt} + par0.preFact(jj)*permute(Sab{tt}(par0.type{jj}(2),par0.type{jj}(3),:,:),[3 4 1 2]);
                    otherwise
                        error('sw:horace:WrongPar','Wrong ''component'' parameter!');
                end
            end
        end
    end

    % normalised volume fractions of the twins
    vol = obj.twin.vol/sum(obj.twin.vol);
    for tt = 1:nTwin
        for ii = 1:size(DSF,1)
            DSF{ii,tt}    = DSF{ii,tt}*vol(tt);
        end
    end

    % add all modes for different twins
    % use only the real part of the dispersion
    omega = real(cell2mat(omega'));
    DSF   = abs(cell2mat(DSF'));

    % normalize intensities
    if param.norm
        % Lande' g-factor
        if spectra.gtensor
            % g-tensor is already included in the spinwave calculation
            g = 1;
        else
            % use simple g=2 here
            g = 2;
        end
        % constant: p = gamma*r0/2
        % neutron magnetic moment constant: M = gamma*gammaN
        gamma = 1.91304272; % 1/s/T
        % classical radius of the electron
        r0 = 2.8179403267e-15; % m
        % cross section constant in mbarn
        p2 = (g*gamma*r0/2)^2*1e28*1e3; % mbarn

        % convert intensity to mbarn/meV units using the energy bin size
        DSF = DSF*p2/param.dE;
    end
    
    for jj = 1:2*nMagExt
        % dispersion in cell
        w{jj}(hklIdxMEM) = omega(jj,:);
        % intensity in cell
        s{jj}(hklIdxMEM) = DSF(jj,:);
    end
end

if param.norm
    fprintf0(obj.fileid,'Intensity is converted to mbarn/meV units.\n');
    if spectra.gtensor
        fprintf0(obj.fileid,'g-tensor was already included in the spin wave calculation.\n');
    else
        fprintf0(obj.fileid,'Isotropic g-tensor of 2 assumed here.\n');
    end
end

end
