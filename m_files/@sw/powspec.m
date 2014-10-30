function spectra = powspec(obj, hklA, varargin)
% calculates powder averaged spin wave spectra
%
% spectra = POWSPEC(obj, hklA, 'Option1', Value1, ...)
%
% Input:
%
% obj       sw class object.
% hklA      Vector containing the Q values in inverse Angstrom where powder
%           spectra will be calculated, dimensions are [1 nQ].
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
%
% 'spectra' is a struct type variable with the following fields:
% swConv    The spectra convoluted with the dispersion. The center
%           of the energy bins are stored in spectra.Evect. Dimensions are
%           [nE nQ].
% hklA      Same Q values as the input hklA [1 nQ]. Evect
%           Contains the input energy transfer values, dimensions are
%           [1 nE].
% param     Contains all the input parameters.
% obj       The copy of the input obj object.
%
% Example:
%
% tri = sw_model('triAF',1);
% E = linspace(0,3,100);
% Q = linspace(0,4,300);
% triSpec = tri.powspec(Q,'Evect',E,'nRand',1e3);
% sw_plotspec(triSpec);
%
% The example calculates the powder spectrum of the triangular lattice
% antiferromagnet (S=1, J=1) between Q = 0 and 3 A^-1 (the lattice
% parameter is 3 Angstrom).
%
% See also SW, SW.SPINWAVE, SW.OPTMAGSTR.
%

% help when executed without argument
if nargin==1
    help sw.powspec
    return
end

fid = obj.fid;

% if function is terminated using Ctrl+C, the original fileid value is
% restored
c = onCleanup(@()obj.fileid(fid));


hklA = hklA(:)';
T0 = obj.single_ion.T;

inpForm.fname  = {'nRand' 'Evect'           'T'   'formfact' 'formfactfun' 'Hermit'};
inpForm.defval = {100     linspace(0,1,100) T0    false      @sw_mff        true    };
inpForm.size   = {[1 1]   [1 -1]            [1 1] [1 -2]     [1 1]          [1 1]   };

param  = sw_readparam(inpForm, varargin{:});

nQ      = length(hklA);
nE      = length(param.Evect);
powSpec = zeros(nE,nQ);

fprintf0(fid,'Calculating powder spectra:\n');

if fid
    sw_status(0,1);
end

for ii = 1:nQ
    rQ  = randn(3,param.nRand);
    Q   = bsxfun(@rdivide,rQ,sqrt(sum(rQ.^2)))*hklA(ii);
    hkl = (Q'*obj.basisvector)'/2/pi;
    
    % no output from spinwave() function
    obj.fileid(0);
    specQ = obj.spinwave(hkl,'fitmode',true,'notwin',true,'Hermit',param.Hermit,...
        'formfact',param.formfact,'formfactfun',param.formfactfun);
    
    % reset output to original value
    obj.fileid(fid);
    specQ = sw_neutron(specQ,'pol',false);
    specQ.obj = obj;
    specQ = sw_egrid(specQ,'Evect',param.Evect,'T',param.T);
    powSpec(:,ii) = sum(specQ.swConv,2)/param.nRand;
    if fid
        sw_status(ii/nQ*100);
    end
end
if fid
    sw_status(100,2);
end

% save different field into spectra
spectra.swConv   = powSpec;
spectra.hklA     = hklA;
spectra.Evect    = param.Evect;
spectra.component = 'Sperp';
spectra.nRand    = param.nRand;
spectra.T        = param.T;
spectra.obj      = copy(obj);
spectra.norm     = false;
spectra.formfact = specQ.formfact;
spectra.incomm   = specQ.incomm;
spectra.helical  = specQ.helical;
spectra.date     = datetime;

% save all input parameters of spinwave into spectra
spectra.param    = specQ.param;

end