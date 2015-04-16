function [w, s] = horace(obj, qh, qk, ql, varargin)
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
% p             Parameters, currently unused.
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
% See also SW, SW.SPINWAVE.
%

if nargin <= 1
    help sw.horace;
    return;
end

inpForm.fname  = {'component'};
inpForm.defval = {'Sperp'    };
inpForm.size   = {[1 -2]     };

warnState = warning('off','sw_readparam:UnreadInput');
param = sw_readparam(inpForm, varargin{:});

% calculate spin wave spectrum
if nargin > 5
    spectra = obj.spinwave([qh(:) qk(:) ql(:)]',varargin{:});
else
    spectra = obj.spinwave([qh(:) qk(:) ql(:)]');
end
warning(warnState);

% calculate Sperp
spectra = sw_neutron(spectra,'pol',false);


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
    Sab   = spectra.Sab;
    Sperp = spectra.Sperp;
    
else
    nTwin = 1;
    omega = {spectra.omega};
    Sab   = {spectra.Sab};
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
                    error('sw:horace:WrongPar','Wrong ''component'' parameter!');
            end
        end
    end
end

% normalised volume fractions of the twins
vol = spectra.obj.twin.vol/sum(spectra.obj.twin.vol);
for tt = 1:nTwin
    for ii = 1:size(DSF,1)
        DSF{ii,tt}    = DSF{ii,tt}*vol(tt);
    end
end

% add all modes for different twins
% use only the real part of the dispersion
omega = real(cell2mat(omega'));
DSF   = abs(cell2mat(DSF'));

% dispersion in cell
w = mat2cell(omega',nHkl,ones(nMode*nTwin,1));
% intensity in cell
s = mat2cell(DSF' ,nHkl,ones(nMode*nTwin,1));

end
