function spectra = sw_instrument(spectra, varargin)
% spectra = SW_INSTRUMENT(spectra, 'option1','value1',...) includes
% instrumental factors (resolution, energy transfer range, etc.) to the
% convoluted spectrum.
%
% Options:
%
% dE            Defines the energy resolution of the instrument. It can be
%               a string, single number or vector.
%                   String  File name, that contains the FWHM energy
%                           resolution values as a function of energy
%                           transfer. The file has to contain two columns,
%                           first is the energy values, the second is the
%                           FWHM resolution at the given energy transfer
%                           value.
%                   Number  Constant FWHM energy resolution as a function
%                           of energy transfer.
%                   Vector  Fitted polynomial of the energy resolution as a
%                           function of energy transfer, see polyfit and
%                           polyval.
% polDeg        Degree of the fitted polynomial to the instrumental
%               resolution data. Default is 5.
% dQ            Momentum transfer resolution of the instrument, FWHM is
%               given in A-1 units, default is 0.
% ThetaMin      Minimum scattering angle in degree, default is 0.
% plot          If the resolution is read from file and plot option is
%               true, tre resolution will be plotted, default is true.
% ki            Momentum of the incident neutrons in A^-1 units.
% Ei            Energy of the incident neutrons in meV.
% formFact      String, contains the name of the magnetic ion in FullProf
%               notation (e.g. Cr^{3+} --> 'MCR3' or 'Cr3'), see sw_mff. It
%               can be also a vector of the 7 coefficients according to the
%               following formula:
%               j0(Qs)> = A*exp(-a*Qs^2) + B*exp(-b*Qs^2) + C*exp(-c*Qs^2) + D,
%               where Qs = Q/(4*pi) and A, a, B, ... are the 7
%               coefficients. Default is 'auto', when the name of the
%               magnetic ion is determined by the first magnetic atom in
%               the crystal.
%
% Output:
%
% spectra       Struct variable, same as input with following additional
%               fields:
%
% formFactor    Magnetic ion name string, or the values of the form factor
%               coefficients, depending on input.
% norm          True, if the spectrum is normalised to mbarn units.
% ki            Incident neutron wave vector as given in the input.
% dE            Energy resolution polynomial as given in the input.
% dQ            FWHM of the momentum resolution.
% swRaw         Original simulated data, withouth the application of the
%               instrumental factors.
%
%
% See also SW_MFF, POLYFIT, POLYVAL.
%

inpForm.fname  = {'dE'  'ki'  'Ei' 'plot' 'polDeg' 'ThetaMin' 'formFact' 'dQ' };
inpForm.defval = {0      0     0     true  5        0          'auto'    0    };
inpForm.size   = {[1 -1] [1 1] [1 1] [1 1] [1 1]    [1 1]      [1 -2]    [1 1]};

param = sw_readparam(inpForm, varargin{:});

if isfield(spectra,'swRaw')
    % take raw convoluted spectra if exists
    spectra.swConv = spectra.swRaw;
else
    % otherwise saw the spectra into swRaw field for future use
    spectra.swRaw = spectra.swConv;
end

% remove non-finite values
spectra.swConv(~isfinite(spectra.swConv)) = 0;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% energy resolution
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if isstring(param.dE)
    
    % load file and read values
    fid = fopen(param.dE);
    res = fscanf(fid,'%f',[2 inf]);
    fclose(fid);
    
    xres = res(1,:);
    yres = res(2,:);
    
    % fit polynom to instrument energy resolution
    polyRes = polyfit(xres,yres,param.polDeg);
    xnew = linspace(min(xres),max(xres),500);
    ynew = polyval(polyRes,xnew);
    
    param.dE = polyRes;
    
    if param.plot
        figure;
        plot(xres,yres,'o-');
        hold all
        plot(xnew,ynew,'r-');
        xlabel('Energy Transfer (meV)');
        ylabel('FWHM energy resolution (meV)');
        title('Polynomial fit the instrumental energy resolution');
        legend('Tabulated resolution',sprintf('Fitted polynomial (degree %d)',param.polDeg));
    end
    
end

if ~any(param.dE)
    param.dE = 1e-8;
else
    spectra.dE = param.dE;
    fprintf('Finite instrumental energy resolution is applied.\n');
end

swConv2 = spectra.swConv * 0;

% convolute gaussian to the convoluted spectra
E = spectra.Evect;

for ii = 1:numel(E)
    % standard deviation of the energy resolution Gaussian
    stdG = polyval(param.dE,E(ii))/2.35482;
    
    % Gaussian with intensity normalised to 1, centered on E(ii)
    fG = exp(-((E-E(ii))/stdG).^2/2);
    fG = fG/sum(fG);
    
    swConv2 = swConv2 + fG' * spectra.swConv(ii,:);
    
end

spectra.swConv = swConv2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Q resolution
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Q = spectra.hklA;

if param.dQ > 0
    
    % standard deviation of the Q resolution Gaussian
    stdG = param.dQ/2.35482;
    
    swConv2 = spectra.swConv * 0;
    for ii = 1:numel(Q)
        % Gaussian with intensity normalised to 1, centered on E(ii)
        fG = exp(-((Q-Q(ii))/stdG).^2/2);
        fG = fG/sum(fG);
        swConv2 = swConv2 + spectra.swConv(:,ii) * fG;
        
    end
    spectra.swConv = swConv2;
    spectra.dQ = param.dQ;
    fprintf('Finite instrumental momentum resolution of %5.3f is applied.\n',param.dQ);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% energy transfer range
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if param.ki <= 0
    param.ki = sw_converter('meV',param.Ei,'k');
end

ki = param.ki;
spectra.ki = ki;

if ki > 0
    cosT = cosd(param.ThetaMin);
    sinT = sind(param.ThetaMin);
    
    Emax = (ki^2-(ki*cosT-sqrt(Q.^2-ki^2*sinT^2)).^2) * sw_converter('k',1,'meV');
    Emin = (ki^2-(ki*cosT+sqrt(Q.^2-ki^2*sinT^2)).^2) * sw_converter('k',1,'meV');
    
    Emax(abs(imag(Emax))>0) = 0;
    Emin(abs(imag(Emin))>0) = 0;
    
    
    %Emax  =  spectra.hklA.*(2*param.ki-spectra.hklA) * sw_converter('k',1,'meV');
    %Emin  = -spectra.hklA.*(2*param.ki+spectra.hklA) * sw_converter('k',1,'meV');
    Elist = repmat(spectra.Evect',[1 size(spectra.swConv,2)]);
    Emin  = repmat(Emin,[size(spectra.swConv,1) 1]);
    Emax  = repmat(Emax,[size(spectra.swConv,1) 1]);
    
    spectra.swConv((Elist<Emin) | (Elist>Emax)) = NaN;
    fprintf('Energy transfer is limited to instrument, using ki=%5.3f A-1.\n',ki);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% magnetic form factor
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if strcmp(param.formFact,'auto')
    obj = spectra.obj;
    if any(obj.unit_cell.S>0)
        mLabel = obj.unit_cell.label(obj.unit_cell.S>0);
        param.formFact = mLabel{1};
    else
        param.formFact = '';
    end
end

% form factor of the given ion
[formFactCalc, formFactCoeff] = sw_mff(param.formFact,spectra.hklA);

if any(formFactCoeff(1:end-1))
    spectra.swConv = bsxfun(@times,spectra.swConv,formFactCalc.^2);
    if isstring(param.formFact)
        fprintf('Magnetic form factor of %s is applied.\n',param.formFact);
    else
        fprintf('Magnetic form factor with given coefficients is applied.\n');
    end
else
    fprintf('No magnetic form factor applied.\n');
end

spectra.formFact = param.formFact;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% normalise spectrum to mbarn/meV
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Lande' g-factor
g  = 2;
% constant: p = gamma*r0/2
% neutron magnetic moment constant: M = gamma*gammaN
gamma = 1.91304272; % 1/s/T
% classical radius of the electron
r0 = 2.8179403267e-15; % m
% cross section constant in mbarn
p2 = (g*gamma*r0/2)^2*1e28*1e3; % mbarn

% convert intensity to mbarn/meV units
%dE = diff(spectra.Evect);
%dE = [dE(1) (dE(2:end)+dE(1:end-1))/2 dE(end)];
spectra.swConv = spectra.swConv*p2;%./repmat(dE',[1 size(spectra.swConv,2)]);
% set 'normalized units' switch on
spectra.norm = true;

fprintf('Intensity is converted to mbarn/meV units.\n');

end
