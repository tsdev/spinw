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
% Evect     Vector, defines the center/edge of the energy bins of the
%           calculated output, dimensions are is [1 nE]. The energy units
%           are defined by the unit.kB property of the sw object. Default
%           value is an edge bin: linspace(0,1.1,101).
% binType   String, determines the type of bin give, possible options:
%               'cbin'    Center bin, the center of each energy bin is given.
%               'ebin'    Edge bin, the edges of each bin is given.
%           Default is 'ebin'.
% T         Temperature to calculate the Bose factor in units
%           depending on the Boltzmann constant. Default is taken from
%           obj.single_ion.T value.
% title     Gives a title string to the simulation that is saved in the
%           output.
% specfun   Function handle of the spectrum calculation function. Default
%           is @spinwave.
% extrap    If true, arbitrary additional parameters are passed over to
%           the spectrum calculation function.
% fibo      If true, instead of random sampling of the unit sphere the
%           Fibonacci numerical integration is implemented as described in:
%           J. Phys. A: Math. Gen. 37 (2004) 11591
%           The number of points on the sphere is given by the largest
%           Fibonacci number below nRand. Default is false.
% useMex        If true, the code will use compiled mex files (if they
%               exist) to speed up the calculation, for details see
%               sw_mex() function. Default is false.
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
% Evect     Energy grid converted to edge bins.
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

title0 = 'Powder LSWT spectrum';

inpForm.fname  = {'nRand' 'Evect'             'T'   'formfact' 'formfactfun'};
inpForm.defval = {100     linspace(0,1.1,101) T0    false      @sw_mff      };
inpForm.size   = {[1 1]   [1 -1]              [1 1] [1 -2]     [1 1]        };

inpForm.fname  = [inpForm.fname  {'hermit' 'gtensor' 'title' 'specfun' }];
inpForm.defval = [inpForm.defval {true     false     title0  @spinwave }];
inpForm.size   = [inpForm.size   {[1 1]    [1 1]     [1 -3]  [1 1]     }];

inpForm.fname  = [inpForm.fname  {'extrap' 'fibo' 'optmem' 'binType' 'useMex'}];
inpForm.defval = [inpForm.defval {false    false  0        'ebin'    false   }];
inpForm.size   = [inpForm.size   {[1 1]    [1 1]  [1 1]    [1 -4]    [1 1]   }];

param  = sw_readparam(inpForm, varargin{:});

% number of bins along energy 
switch param.binType
    case 'cbin'
        nE      = numel(param.Evect);
    case 'ebin'
        nE      = numel(param.Evect) - 1;
end

nQ      = numel(hklA);
powSpec = zeros(nE,nQ);

fprintf0(fid,'Calculating powder spectra:\n');

if fid
    sw_status(0,1);
end

if param.fibo
    % apply the Fibonacci numerical integration on a sphere
    % according to J. Phys. A: Math. Gen. 37 (2004) 11591
    % create QF points on the unit sphere
    
    [F,F1] = sw_fibo(param.nRand);
    param.nRand = F;
    
    QF = zeros(3,F);
    
    j = 0:(F-1);
    QF(3,:) = j*2/F-1;
    
    theta = asin(QF(3,:));
    phi   = 2*pi*F1/F*j;
    
    QF(1,:) = cos(theta).*sin(phi);
    QF(2,:) = cos(theta).*cos(phi);
    
end

for ii = 1:nQ
    if param.fibo
        Q = QF*hklA(ii);
    else
        rQ  = randn(3,param.nRand);
        Q   = bsxfun(@rdivide,rQ,sqrt(sum(rQ.^2)))*hklA(ii);
    end
    hkl = (Q'*obj.basisvector)'/2/pi;
    
    % no output from spinwave() function
    obj.fileid(0);
    if param.extrap
        % allow arbitrary additional parameters to pass to the spectral
        % calculation function
        warnState = warning('off','sw_readparam:UnreadInput');
        specQ = param.specfun(obj,hkl,varargin{:});
        warning(warnState);
        
    else
        specQ = param.specfun(obj,hkl,'fitmode',true,'notwin',true,...
            'Hermit',param.hermit,'formfact',param.formfact,...
            'formfactfun',param.formfactfun,'gtensor',param.gtensor,...
            'optmem',param.optmem,'useMex',param.useMex);
    end
    
    % reset output to original value
    obj.fileid(fid);
    specQ = sw_neutron(specQ,'pol',false);
    specQ.obj = obj;
    % use edge grid by default
    specQ = sw_egrid(specQ,'Evect',param.Evect,'T',param.T,'binType',param.binType);
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
spectra.Evect    = specQ.Evect;
spectra.component = 'Sperp';
spectra.nRand    = param.nRand;
spectra.T        = param.T;
spectra.obj      = copy(obj);
spectra.norm     = false;
spectra.formfact = specQ.formfact;
spectra.gtensor  = specQ.gtensor;
spectra.incomm   = specQ.incomm;
spectra.helical  = specQ.helical;
spectra.date     = datestr(now);
spectra.title    = param.title;

% save all input parameters of spinwave into spectra
spectra.param    = specQ.param;

end