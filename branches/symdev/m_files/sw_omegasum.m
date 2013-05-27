function spectra = sw_omegasum(spectra, varargin)
% spec = sw_omegasum(spec, 'Option1', Value1, ...) removes the degenerate
% modes from the dispersion stored in spectra.omega and sorts omega
% according to the energy.
%
% The degenerate dispersion energies are substituted with NaN values. Be
% carefull, after this function sw_conv won't work properly on spectra.
%
% Options:
%
% tol       Tolerance, within two energies are considered equal. Default
%           value is 1e-5.
% zeroint   The minimum intensity value, below the mode is dropped. Default
%           value is zero (no modes are dropped due to weak intensity).
%

inpForm.fname  = {'tol' 'zeroint'};
inpForm.defval = {1e-5  0        };
inpForm.size   = {[1 1] [1 1]    };

param = sw_readparam(inpForm, varargin{:});

omega    = abs(real(spectra.omega));
omega(isnan(omega)) = 0;
omegaCol = omega*NaN;

nQ = size(omega,2);

if isfield(spectra,'swInt') && param.zeroint > 0
    omega(abs(real(spectra.swInt)) < param.zeroint) = 0;
end

for ii = 1:nQ
    oTemp = sw_uniquetol(omega(:,ii)',param.tol);
    oTemp(oTemp == 0) = NaN;
    omegaCol(1:numel(oTemp),ii) = sort(oTemp);
end

spectra.omega = omegaCol;

end