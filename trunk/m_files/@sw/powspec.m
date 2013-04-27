function spectra = powspec(obj, hklA, varargin)
% calculates powder averaged spin wave spectra
%
% spectra = POWSPEC(obj, hklA, 'Option1', Value1, ...)
%
% hklA          Q values in reciproc angstrom, where powder spectra is
%               to be calculated, dimensions are[1 nQ].
%
% Options:
%
% swfunc        Function for the spin wave calculation. Default is
%               <a href="matlab: doc sw.spinwave">@spinwave</a>.
% nRand         Number of random orientations per Q value, default is 100.
% Evect         Vector, defined the energy transfer values for the
%               convoluted output in units of meV, dimensions are [1 nE].
%               Default is linspace(0,1,100).
% T             Temperature to calculate the Bose factor in units
%               depending on the Boltzmann constant. Default is zero.
% formfact      TODO!!!!!!!!!!!!!!!!!!!!!!!!
%               Whether to include the magnetic form factor in the swConv
%               convoluted spectra. If true the form factor based on the
%               name of the magnetic atoms are used to read form factor
%               data from formfactor.dat file.
%
% Output:
% spectra is struct type with the following fields:
%
% swConv        The spectra convoluted with the dispersion. The center
%               of the energy bins are stored in spectra.Evect. Dimensions
%               are [nE nQ].
% hklA          Same Q values as the input hklA [1 nQ]. Evect
%               Contains the input energy transfer values, dimensions are
%               [1 nE].
% param         Contains all the input parameters.
% obj           The input sw object.
%
% See also SW, SW.SWINC, SW.SPINWAVE, SW.OPTMAGSTR.
%

% help when executed without argument
if nargin==1
    help sw.powspec
    return
end

hklA = hklA(:)';

inpForm.fname  = {'swfunc'  'nRand' 'Evect'           'T'  };
inpForm.defval = {@spinwave 100     linspace(0,1,100) 0    };
inpForm.size   = {[1 1]     [1 1]   [1 -1]            [1 1]};

param  = sw_readparam(inpForm, varargin{:});

nQ      = length(hklA);
nE      = length(param.Evect);
powSpec = zeros(nE,nQ);

fprintf('Calculating powder spectra:\n');
sw_status(0,1);
for ii = 1:nQ
    sw_status(ii/nQ*100);
    rQ  = randn(3,param.nRand);
    Q   = bsxfun(@rdivide,rQ,sqrt(sum(rQ.^2)))*hklA(ii);
    hkl = (Q'*obj.basisvector)'/2/pi;
    
    specQ = param.swfunc(obj,hkl,'fitmode',false,'notwin',true);
    specQ = sw_neutron(specQ,'pol',false);
    specQ = sw_conv(specQ,'Evect',param.Evect,'T',param.T);
    powSpec(:,ii) = sum(specQ.swConv,2)/param.nRand;
end
sw_status(100,2);

spectra.swConv   = powSpec;
spectra.hklA     = hklA;
spectra.Evect    = param.Evect;
spectra.swfunc   = param.swfunc;
spectra.convmode = 'Sperp';
spectra.nRand    = param.nRand;
spectra.T        = param.T;
spectra.obj      = copy(obj);

end