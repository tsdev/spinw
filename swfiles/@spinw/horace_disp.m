function [w, s] = horace(obj, qh, qk, ql, varargin)
% calculates spin wave dispersion/correlation functions to be called from Horace
%
% [w, s] = HORACE(obj, qh, qk, ql, 'Option1', Value1, ...)
%
% The function produces spin wave dispersion and intensity for Horace
% (<a href=http://horace.isis.rl.ac.uk>http://horace.isis.rl.ac.uk</a>).
%
% Input:
%
% obj           Input spinw object.
% qh, qk, ql    Reciprocal lattice components in reciprocal lattice units.
%
% Options:
%
% component Selects the previously calculated intensity component to be
%           convoluted. The possible options are:
%               'Sperp' convolutes the magnetic neutron scattering
%                       intensity (<Sperp * Sperp> expectation value).
%                       Default.
%               'Sab'   convolutes the selected components of the spin-spin
%                       correlation function. Letter a and b can be 'x',
%                       'y' or 'z'. For example: 'Sxx' will convolute the
%                       xx component of the correlation function with the
%                       dispersion. xyz is the standard coordinate system,
%                       see online documentation of SpinW.
%           Any linear combination of the above are allowed, for example:
%           'Sxx+2*Syy' convolutes the linear combination of the xx
%           component of the spin-spin correlation function and the yy
%           component.
% norm      If true the spin wave intensity is normalized to mbarn/meV/(unit
%           cell) units. Default is false.
% dE        Energy bin size, for intensity normalization. Use 1 for no
%           division by dE in the intensity.
% param     Input parameters (can be used also within Tobyfit). Additional
%           options ('mat','selector') might be necessary, for details see
%           sw.matparser function. All extra parameters of sw.horace
%           function will be forwarded to the sw.matparser function before
%           calculating the spin wave spectrum (or any user written parser
%           function). For user written functions defined with the
%           following header:
%               func(obj,param)
%           the value of the param option will be forwarded. For user
%           functions with variable number of arguments, all input options
%           of sw.horace will be forwarded. In this case it is recommended
%           to use sw_readparam() function to handle the variable number
%           arguments within func().
% parfunc   Parser function of the 'param' input. Default is
%           @sw.matparser which can be used directly by Tobyfit. For user
%           defined functions the minimum header has to be:
%               func(obj,param)
%           where obj is an spinw type object, param is the parameter
%           values forwarded from spinw.horace directly.
% func      User function that will be called after the parameters set on
%           the SpinW object. It can be used to optimize magnetic
%           structure for the new parameters, etc. The input should be a
%           function handle of a function with a header:
%               fun(obj)
%
% useFast   whether to use the SPINW.SPINWAVEFAST method or not. This method
%           is similar to SPINW.SPINWAVE but calculates only the unpolarised
%           neutron cross-section and ignores all negative energy branches
%           as well as using other shortcuts. In general it should produce
%           the same spectra as SPINW.SPINWAVE, with some rounding errors,
%           but can be 2-3 times faster and uses less memory. 
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
% See also SPINW, SPINW.SPINWAVE, SPINW.MATPARSER, SW_READPARAM.
%

if nargin <= 1
    help spinw.horace_disp
    return
end

inpForm.fname  = {'component' 'norm' 'dE'  'parfunc'      'param' 'func' 'useFast'};
inpForm.defval = {'Sperp'     false  0     @obj.matparser []      []     true};
inpForm.size   = {[1 -1]      [1 1]  [1 1] [1 1]          [1 -2]  [1 1]  [1 1]};
inpForm.soft   = {false       false  false false          true    true   false};

warnState = warning('off','sw_readparam:UnreadInput');
% To handle matlab fitting syntax, which may set the first argument to the
% fittable parameter (without keyword).
if ~ischar(varargin{1})
    param = sw_readparam(inpForm, varargin{2:end});
    if isempty(param.param)
        varargin = {varargin{:} 'param' varargin{1}};
    end
    param.param = varargin{1};  % Always override keyword specified parameters
    varargin(1) = [];
else
    param = sw_readparam(inpForm, varargin{:});
end

if ~isempty(param.param)
    % change matrix values for Horace data fitting
    if nargin(param.parfunc) < 0
        param.parfunc(varargin{:});
    elseif nargin(param.parfunc) == 2
        param.parfunc(param.param);
    else
        error('spinw:horace:WrongInput','User defined function with incompatible header!');
    end
end
%}

if param.norm && param.dE == 0
    error('spinw:horace:WrongInput',['To convert spin wave intensity to mbarn/meV/cell/sr'...
        ' units, give the energy bin step.'])
end

% call user defined function on the spinw object
if ~isempty(param.func)
    param.func(obj);
end

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

needSab = false;
for ii = 1:numel(parsed)
    par0 = parsed{ii};
    for jj = 1:length(par0.type)
        if par0.type{jj}(1) ~= 1
            needSab = true;
            break;
        end
    end
end

% calculate spin wave spectrum
if nargin > 5
    % include the fitmode option to speed up calculation
    if numel(varargin) == 1
        varargin{1}.fitmode = 2;
        if needSab || ~param.useFast
            spectra = obj.spinwave([qh(:) qk(:) ql(:)]',varargin{1});
        else
            spectra = obj.spinwavefast([qh(:) qk(:) ql(:)]',varargin{1});
        end
    else
        if needSab || ~param.useFast
            spectra = obj.spinwave([qh(:) qk(:) ql(:)]',varargin{:},'fitmode',true);
        else
            spectra = obj.spinwavefast([qh(:) qk(:) ql(:)]',varargin{:},'fitmode',true);
        end
    end
else
    if needSab || ~param.useFast
        spectra = obj.spinwave([qh(:) qk(:) ql(:)]','fitmode',true);
    else
        spectra = obj.spinwavefast([qh(:) qk(:) ql(:)]','fitmode',true);
    end
end
warning(warnState);

% calculate Sperp
if needSab || ~param.useFast
    spectra = sw_neutron(spectra,varargin{:},'pol',false);
end


% pack all cross section into a cell for easier looping
if iscell(spectra.omega)
    nTwin = numel(spectra.omega);
    omega = spectra.omega;
    if needSab
        Sab   = spectra.Sab;
    end
    Sperp = spectra.Sperp;
    
else
    nTwin = 1;
    omega = {spectra.omega};
    if needSab
        Sab   = {spectra.Sab};
    end
    Sperp = {spectra.Sperp};
end

% extract the requested cross section
nMode = size(omega{1},1);
nHkl  = size(omega{1},2);

% DSF stores the intensity that is going to be convoluted
DSF = cell(nConv,nTwin);
% select value from Sperp or Sab matrices
for tt = 1:nTwin
    for ii = 1:numel(parsed)
        par0 = parsed{ii};
        DSF{ii,tt} = zeros(nMode, nHkl);
        for jj = 1:length(par0.type)
            switch par0.type{jj}(1)
                case 1
                    DSF{ii,tt} = DSF{ii,tt} + par0.preFact(jj)*Sperp{tt};
                case 2
                    DSF{ii,tt} = DSF{ii,tt} + par0.preFact(jj)*permute(Sab{tt}(par0.type{jj}(2),par0.type{jj}(3),:,:),[3 4 1 2]);
                otherwise
                    error('spinw:horace:WrongPar','Wrong ''component'' parameter!');
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
    
    fprintf0(obj.fileid,'Intensity is converted to mbarn/meV units.\n');
    if spectra.gtensor
        fprintf0(obj.fileid,'g-tensor was already included in the spin wave calculation.\n');
    else
        fprintf0(obj.fileid,'Isotropic g-tensor of 2 assumed here.\n');
    end
end

% dispersion in cell
w = mat2cell(omega',nHkl,ones(nMode*nTwin,1));
% intensity in cell
s = mat2cell(DSF' ,nHkl,ones(nMode*nTwin,1));

end
