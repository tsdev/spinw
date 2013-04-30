function [omega, swConv] = sw_idata(spectrum, varargin)
% [omega, swConv] = SW_IDATA(spectrum, 'option1', value1 ...) creates iData object
% (<a href=http://ifit.mccode.org>http://ifit.mccode.org</a>) and convolutes the spectra with a fixed
% instrumental resolution, assuming the energy and Q axis are linear.
%
% Input:
%
% spectrum      Calculated spin wave spectrum, struct type.
%
% Options:
%
% fwhmE         Full width half maximum of the Gaussian energy
%               resolution. Works properly only for linear energy axis.
% fwhmQ         Full width half maximum of the Gaussian momentum
%               transfer resolution in A^-1 units. Works properly only
%               for linear scans in reciprocal space. Be carefull for
%               example for scans like [1, QK, 0] where the equal QK
%               steps give unequal steps in the A^-1 reciprocal space.
% nInterp       Number of axis subdivision before convolution, equal
%               for Q and E.
%


inpForm.fname  = {'fwhmE' 'fwhmQ' 'nInterp'};
inpForm.defval = {0.1     0.01    1        };
inpForm.size   = {[1 1]   [1 1]   [1 1]    };

param = sw_readparam(inpForm, varargin{:});


iDataInstalled = exist('iData') && isa(iData,'iData'); %#ok<EXIST>

if ~iDataInstalled
    error('sw:sw_idata:NoIdataClass','iFit is not installed, install it from: http://ifit.mccode.org')
end

% Labels for the x-axis
if isfield(spectrum,'hkl')
    [xLabel, xAxis] = sw_label(spectrum);
else
    % Powder mode
    xLabel = 'Momentum Transfer (A^-1)';
    xAxis  = spectrum.hklA;
end

% Save omega values as iData array
if isfield(spectrum,'omega')
    nMagExt = size(spectrum.omega,1)/2;
    
    % Save the absolute real part of the dispersion
    omega = iData(1:2*nMagExt,xAxis,abs(real(spectrum.omega)'));
    
    omega.Error = omega.Signal*0;
    
    label(omega,0,'Energy transfer (meV)');
    xlabel(omega,'Mode index');
    ylabel(omega,xLabel);
    
    omega.Title = 'Spin wave spectra';
    omega = dog(2,omega);
else
    omega = [];
end

% Save the convoluted dispersion
if isfield(spectrum,'swConv')
    swConv = iData(xAxis,spectrum.Evect,abs(real(spectrum.swConv)));
    swConv.Error = swConv.Signal*0;
    label(swConv,0,'Intensity (arb. units)');
    xlabel(swConv, xLabel);
    ylabel(swConv,'Energy transfer (meV)');
    swConv.Title = 'Convoluted spin wave spectra';
else
    swConv = iData;
end


% Convolute data with resolution function assuming equal energy steps
stepE = sum(spectrum.Evect(end)-spectrum.Evect(1))/(length(spectrum.Evect)-1)/param.nInterp;
% Average Q steps
absQ  = sqrt(sum(spectrum.hklA.^2,1));
stepQ = sum(abs(absQ(2:end)-absQ(1:(end-1))))/(length(absQ)-1)/param.nInterp;

% Create the 2D Gaussian convolution core
nGx = ceil(param.fwhmE*3/stepE);
nGy = ceil(param.fwhmQ*3/stepQ);
xG  = (-nGx:nGx)*stepE;
yG  = (-nGy:nGy)*stepQ;

[xxG,yyG] = ndgrid(xG,yG);

% Parameters for the 2D Gaussian
pG  = [1 0 0 param.fwhmE/2.35482 param.fwhmQ/2.35482 0 0];

swConv = conv(interp(swConv,param.nInterp),gauss2d(pG,xxG,yyG),'same pad normalize');

end