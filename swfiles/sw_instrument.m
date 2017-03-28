function spectra = sw_instrument(spectra, varargin)
% includes instrumental factors into the calculated spectrum
%
% spectra = SW_INSTRUMENT(spectra, 'option1','value1',...)
%
% It includes instrumental factors (resolution, energy transfer range,
% etc.) to the convoluted spectrum.
%
% Options:
%
% dE            Defines the FWHM energy resolution of the instrument. It
%               can be a string, single number, vector of function hangle:
%                 String    File name, that contains the FWHM energy
%                           resolution values as a function of energy
%                           transfer. The file has to contain two columns,
%                           first is the energy values, the second is the
%                           FWHM resolution at the given energy transfer
%                           value, see sw_res() function for details.
%                 Number    Constant FWHM energy resolution as a function
%                           of energy transfer.
%                 Matrix    Dimensions of Nx2, first column contains the
%                           energy transfer values, second column contains
%                           the FWHM resolution values. These values will
%                           be fitted using a polynomial with a fixed
%                           degree, see sw_res() function for details.
%                 Function  Function handle of a resolution function
%                           with the following header:
%                               E_FWHM = res_fun(E)
%                           where E_FWHM is the FWHM energy resolution and
%                           E is the energy transfer value.
% func          Shape of the energy resolution function, for details see
%               the help of sw_resconv.
% polDeg        Degree of the fitted polynomial to the instrumental
%               resolution data. Default is 5.
% dQ            Momentum transfer resolution of the instrument, FWHM is
%               given in A-1 units, default is 0.
% ThetaMin      Minimum scattering angle in degree, default is 0.
% plot          If the resolution is read from file and plot option is
%               true, tre resolution will be plotted, default is true.
%
% Fixed incident neutron energy:
% ki            Momentum of the incident neutrons in A^-1 units.
% Ei            Energy of the incident neutrons in meV.
%
% Fixed final neutron energy:
% kf            Final momentum of the neutrons in A^-1 units.
% Ef            Final neutron energy in meV.
%
% norm          If true, the data is normalized to mbarn units. Default is
%               false. If no g-tensor is included in the spin wave
%               calculation, g-tensor = 2 is assumed here.
% useRaw        If false, the already modified spectra.swConv field is
%               modified further instead of the original powder spectrum
%               stored in spectra.swRaw. Default is true.
%
% Output:
%
% spectra       Struct variable, same as input with following additional
%               fields:
%
% norm          True, if the spectrum is normalised to mbarn units.
% ki            Incident neutron wave vector as given in the input.
% dE            Energy resolution polynomial as given in the input.
% dQ            FWHM of the momentum resolution.
% swRaw         Original simulated data, withouth the application of the
%               instrumental factors.
%
%
% See also POLYFIT, POLYVAL, SW_RES, SW_RESCONV.
%

if nargin == 0
    help sw_instrument
    return
end

func0 = @swfunc.gaussfwhm;

inpForm.fname  = {'dE'    'ki'  'Ei'  'kf'  'Ef'  'plot' 'polDeg' 'ThetaMin'};
inpForm.defval = {0       0     0     0     0     false   5        0        };
inpForm.size   = {[-1 -2] [1 1] [1 1] [1 1] [1 1] [1 1]  [1 1]    [1 1]     };

inpForm.fname  = [inpForm.fname  {'formFact' 'dQ'  'norm' 'useRaw' 'func'}];
inpForm.defval = [inpForm.defval { 'auto'    0     false   true    func0 }];
inpForm.size   = [inpForm.size   { [1 -2]    [1 1] [1 1]   [1 1]   [1 1] }];

param = sw_readparam(inpForm, varargin{:});

% Print output
fid0 = spectra.obj.fileid;

if isfield(spectra,'swRaw')
    % take raw convoluted spectra if exists and request
    if param.useRaw
        spectra.swConv = spectra.swRaw;
    end
else
    % otherwise saw the spectra into swRaw field for future use
    spectra.swRaw = spectra.swConv;
end

% remove non-finite values
if ~iscell(spectra.swConv)
    spectra.swConv = {spectra.swConv};
end

% number of plots
nPlot = numel(spectra.swConv);

for ii = 1:nPlot
    spectra.swConv{ii}(~isfinite(spectra.swConv{ii})) = 0;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% energy resolution
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% determine whether to include resolution
calcres = true;
try %#ok<TRYNC>
    if ~any(param.dE)
        calcres = false;
    end
end

if isa(param.dE,'function_handle')
    % do nothing
elseif numel(param.dE)>1
    
    % determine the energy resolution curve from a file or the given matrix
    param.dE = sw_res(param.dE,param.polDeg,param.plot);
    
end

% center bin positions
cEvect = (spectra.Evect(1:(end-1))+spectra.Evect(2:end))/2;

% include resolution
if calcres
    for jj = 1:nPlot
        %spectra.swConv{jj} = sw_resconv(spectra.swConv{jj},cEvect',param.dE);
        spectra.swConv{jj} = sw_resconv(spectra.swConv{jj},cEvect',param.dE,param.func);
    end
    
    fprintf0(fid0,'Finite instrumental energy resolution is applied.\n');
end
spectra.dE = param.dE;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Q resolution
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Q = spectra.hklA;

if size(Q,1) > 1
    %Q = sqrt(sum(bsxfun(@minus,Q,Q(:,1)).^2,1));
    Q = sqrt(sum(Q.^2,1));
end

if param.dQ > 0
    
    % standard deviation of the Q resolution Gaussian
    stdG = param.dQ/2.35482;
    
    for jj = 1:nPlot
        swConv = spectra.swConv{jj};
        swConvTemp = swConv * 0;
        for ii = 1:numel(Q)
            % Gaussian with intensity normalised to 1, centered on E(ii)
            fG = exp(-((Q-Q(ii))/stdG).^2/2);
            fG = fG/sum(fG);
            swConvTemp = swConvTemp + swConv(:,ii) * fG;
            
        end
        spectra.swConv{jj} = swConvTemp;
    end
    
    fprintf0(fid0,'Finite instrumental momentum resolution of %5.3f A-1 is applied.\n',param.dQ);
end

spectra.dQ = param.dQ;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% energy transfer range
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if (param.ki>0) || (param.Ei>0)
    FX = 1;
    param.k = param.ki;
    param.E = param.Ei;
elseif (param.kf>0) || (param.Ef>0)
    FX = 2;
    param.k = param.kf;
    param.E = param.Ef;
else
    FX = 0;
end

if FX>0 && param.k == 0
    param.k = sw_converter(param.E,'meV','k');
end

% save ki/kf value into the output
switch FX
    case 1
        spectra.ki = param.k;
        kstr = 'ki';
    case 2
        spectra.kf = param.k;
        kstr = 'kf';
end

if FX > 0
    k0 = param.k;
    cosT = cosd(param.ThetaMin);
    sinT = sind(param.ThetaMin);
    
    
    for jj = 1:nPlot
        switch FX
            case 1
                % fix ki
                Emax = (k0^2-(k0*cosT-sqrt(Q.^2-k0^2*sinT^2)).^2) * sw_converter(1,'k','meV');
                Emin = (k0^2-(k0*cosT+sqrt(Q.^2-k0^2*sinT^2)).^2) * sw_converter(1,'k','meV');
            case 2
                % fix kf
                Emax = -(k0^2-(k0*cosT+sqrt(Q.^2-k0^2*sinT^2)).^2) * sw_converter(1,'k','meV');
                Emin = -(k0^2-(k0*cosT-sqrt(Q.^2-k0^2*sinT^2)).^2) * sw_converter(1,'k','meV');
        end
        %Emax = (ki^2-(ki*cosT-sqrt(Q.^2-ki^2*sinT^2)).^2) * sw_converter(1,'k','meV');
        %Emin = (ki^2-(ki*cosT+sqrt(Q.^2-ki^2*sinT^2)).^2) * sw_converter(1,'k','meV');
        
        Emax(abs(imag(Emax))>0) = 0;
        Emin(abs(imag(Emin))>0) = 0;
        
        Elist = repmat(cEvect',[1 size(spectra.swConv{jj},2)]);
        Emin  = repmat(Emin,[size(spectra.swConv{jj},1) 1]);
        Emax  = repmat(Emax,[size(spectra.swConv{jj},1) 1]);
        
        
        idx = (Elist<Emin) | (Elist>Emax);
        swConv = spectra.swConv{jj};
        swConv(idx) = NaN;
        spectra.swConv{jj} = swConv;
    end
    fprintf0(fid0,'Energy transfer is limited to instrument, using %s=%5.3f A-1.\n',kstr,k0);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% normalise spectrum to mbarn/meV/f.u. or mbarn/meV/cell
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~isfield(spectra,'gtensor')
    spectra.gtensor = false;
end

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
    dE = diff(spectra.Evect);
    % the new Evect is in edge bin mode
    %dE = [dE(1) (dE(2:end)+dE(1:end-1))/2 dE(end)];
    
    for jj = 1:nPlot
        spectra.swConv{jj} = spectra.swConv{jj}*p2./repmat(dE',[1 size(spectra.swConv{jj},2)]);
    end
    
    % set 'normalized units' switch on
    spectra.norm = true;
    if spectra.obj.unit.nformula > 0
        fprintf0(fid0,'Intensity is converted to mbarn/meV/f.u. units.\n');
    else
        fprintf0(fid0,'Intensity is converted to mbarn/meV/cell units.\n');
    end
    if spectra.gtensor
        fprintf0(fid0,'g-tensor was already included in the spin wave calculation.\n');
    else
        fprintf0(fid0,'Isotropic g-tensor of 2 assumed here.\n');
    end
    
else
    spectra.norm = false;
end

if nPlot == 1
    spectra.swConv = spectra.swConv{1};
end

end
