function [w, s] = horace(obj, qh, qk, ql, varargin)
% dispersion/correltion function calculator, can be called from Horace
%
% [w, s] = HORACE(obj, qh, qk, ql, p) function to produce
% spin wave dispersion and intensity for Horace (<a href=http://horace.isis.rl.ac.uk>http://horace.isis.rl.ac.uk</a>).
%
% Input:
%
% ojb           Input sw object.
% qh, qk, ql    Reciprocal lattice components in reciprocal lattice units.
% p             Prameters, not used.
%
% Options:
%
% convmode  Selects the previously calculated intensity component to be
%           convoluted. The possible options are:
%               'Sperp' convolutes the magnetic neutron scattering
%                       intensity (<Sperp * Sperp> expectation value).
%                       Default.
%               'Sab'   convolutes the selected components of the spin-spin
%                       correlation function. Letter a and b can be 'x',
%                       'y' or 'z'. For example: 'Sxx' will convolute the
%                       xx component of the correlation function with the
%                       dispersion. xyz is the standard coordinate system,
%                       see online documentation of sw.
%           Any linear combination of the above are allowed, for example:
%           'Sxx+2*Syy' convolutes the linear combination of the xx
%           component of the spin-spin correlation function and the yy
%           component.
%
% Example:
%
% horace_on;
% d3dobj = d3d(cryst.abc,[0 1 0 0],[0,0.01,1],[0 0 1 0],[0,0.01,1],[0 0 0 1],[0,0.1,10]);
% d3dobj = disp2sqw_eval(d3dobj,@cryst.horace,{'convmode','Sperp'},0.1);
% plot(d3dobj);
%
% This example creates a d3d object, a square in (h,k,0) plane and in
% energy between 0 and 10 meV. Then calculates the neutron scattering
% intensity of the spin wave in this volume and plots it using sliceomatic.
%
% See also SW, SW.SPINWAVE.
%

if nargin <= 1
    help sw.horace;
    return;
end

inpForm.fname  = {'convmode'};
inpForm.defval = {'Sperp'   };
inpForm.size   = {[1 -2]    };

param = sw_readparam(inpForm, varargin{:});

% calculate spin wave spectrum
if nargin > 5
    spectra = obj.spinwave([qh(:) qk(:) ql(:)]',varargin{:},'showWarn',false);
else
    spectra = obj.spinwave([qh(:) qk(:) ql(:)]');
end

% calculate Sperp
spectra = sw_neutron(spectra,'pol',false);


% parse the convmode string
if iscell(param.convmode)
    nConv = numel(param.convmode);
    parsed = cell(1,nConv);
    for ii = 1:numel(param.convmode)
        parsed{ii} = sw_parstr(param.convmode{ii});
    end
elseif isstruct(param.convmode)
    nConv  = 1;
    parsed = {param.convmode};
    param.convmode = {parsed{1}.string};
else
    nConv = 1;
    parsed = {sw_parstr(param.convmode)};
    param.convmode = {param.convmode};
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
                    error('sw:horace:WrongPar','Wrong convmode parameter!');
            end
        end
    end
end


% add all modes for different twins
omega = cell2mat(omega');
DSF = abs(cell2mat(DSF'));

% dispersion in cell
w = mat2cell(omega',nHkl,ones(nMode,1));
% intensity in cell
s = mat2cell(DSF' ,nHkl,ones(nMode,1));

end