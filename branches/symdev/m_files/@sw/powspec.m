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
% nRand     Number of random orientations per Q value, default is 100.
% Evect     Vector, defined the energy transfer values for the
%           convoluted output in units of meV, dimensions are [1 nE].
%           Default is linspace(0,1,100).
% T         Temperature to calculate the Bose factor in units
%           depending on the Boltzmann constant. Default is taken from
%           obj.single_ion.T value.
%
% Output:
% spectra is struct type with the following fields:
%
% swConv    The spectra convoluted with the dispersion. The center
%           of the energy bins are stored in spectra.Evect. Dimensions are
%           [nE nQ].
% hklA      Same Q values as the input hklA [1 nQ]. Evect
%           Contains the input energy transfer values, dimensions are
%           [1 nE].
% param     Contains all the input parameters.
% obj       The input sw object.
%
% See also SW, SW.SPINWAVE, SW.OPTMAGSTR.
%

% help when executed without argument
if nargin==1
    help sw.powspec
    return
end

hklA = hklA(:)';
T0 = obj.single_ion.T;

inpForm.fname  = {'nRand' 'Evect'           'T'   'formfact' 'Hermit'};
inpForm.defval = {100     linspace(0,1,100) T0    false      true    };
inpForm.size   = {[1 1]   [1 -1]            [1 1] [1 -2]     [1 1]   };

param  = sw_readparam(inpForm, varargin{:});

nQ      = length(hklA);
nE      = length(param.Evect);
powSpec = zeros(nE,nQ);

fprintf('Calculating powder spectra:\n');
sw_status(0,1);
for ii = 1:nQ
    rQ  = randn(3,param.nRand);
    Q   = bsxfun(@rdivide,rQ,sqrt(sum(rQ.^2)))*hklA(ii);
    hkl = (Q'*obj.basisvector)'/2/pi;
    
    specQ = obj.spinwave(hkl,'fitmode',true,'notwin',true,'fid',0,'Hermit',param.Hermit);
    specQ = sw_neutron(specQ,'pol',false);
    specQ.obj = obj;
    specQ = sw_egrid(specQ,'Evect',param.Evect,'T',param.T);
    powSpec(:,ii) = sum(specQ.swConv,2)/param.nRand;
    sw_status(ii/nQ*100);
end
sw_status(100,2);

spectra.swConv   = powSpec;
spectra.hklA     = hklA;
spectra.Evect    = param.Evect;
spectra.component = 'Sperp';
spectra.nRand    = param.nRand;
spectra.T        = param.T;
spectra.obj      = copy(obj);
spectra.norm     = false;

end