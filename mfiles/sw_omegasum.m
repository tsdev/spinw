function spectra = sw_omegasum(spectra, varargin)
% removes degenerate and ghost magnon modes from spectrum
%
% spec = sw_omegasum(spec, 'Option1', Value1, ...)
%
% It removes the degenerate modes from the dispersion stored in
% spectra.omega and sorts omega according to the energy.
%
% The degenerate dispersion energies are substituted with NaN values. Be
% carefull, after this function sw_egrid() won't work properly on spectra.
% It doesn't work for spectra with multiple twins.
%
% Options:
%
% tol       Tolerance, within two energies are considered equal. Default
%           value is 1e-5.
% zeroint   The minimum intensity value, below the mode is dropped. Default
%           value is zero (no modes are dropped due to weak intensity).
%
% See also SW.SPINWAVE, SW_EGRID.
%

if iscell(spectra.omega)
    error('sw_omegasum:NoTwin','The sw_omegasum() function doesn''t work for spectra calculated for multiple twins!');
end

inpForm.fname  = {'tol' 'zeroint'};
inpForm.defval = {1e-5  0        };
inpForm.size   = {[1 1] [1 1]    };

param = sw_readparam(inpForm, varargin{:});

tol = param.tol;
omega    = real(spectra.omega);
omega(isnan(omega)) = 0;
omega(omega<0) = 0;
omegaCol = omega*NaN;
intCol   = omegaCol;

nQ = size(omega,2);

if isfield(spectra,'swInt') && param.zeroint > 0
    omega(abs(real(spectra.swInt)) < param.zeroint) = 0;
end

for ii = 1:nQ
    oTemp = sw_uniquetol(omega(:,ii)',param.tol);
    oTemp(oTemp == 0) = NaN;
    omegaCol(1:numel(oTemp),ii) = sort(oTemp);
    oTemp(isnan(oTemp)) = [];
    
    if isfield(spectra,'swInt') && ~isempty(oTemp)
        for jj = 1:numel(oTemp)
            intCol(jj,ii) = sum(spectra.swInt(abs(omegaCol(jj,ii)-omega(:,ii))<tol,ii));
        end
    end
    
end

spectra.omega = omegaCol;
spectra.swInt = intCol;
end