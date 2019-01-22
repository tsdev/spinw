function weight = horace_sqw(obj, qh, qk, ql, en, pars, varargin)
% Calculate spectral weight from a spinW model for Horace. Uses disp2sqw
% as the back-end function to calculate the convolution.
%
%   >> weight = swobj.horace_sqw(qh,qk,ql,en,pars,swobj,pars,kwpars)
%
% Input:
% ------
%   qh,qk,ql,en Arrays containing points at which to evaluate sqw from the
%               broadened dispersion
%
%   pars        Arguments needed by the function.
%               - pars = [model_pars scale_factor resolution_pars]
%               - Should be a vector of parameters
%               - The first N parameters relate to the spin wave dispersion
%                 and correspond to spinW matrices in the order defined by
%                 the 'mat' option [N=numel(mat)]
%               - The next M parameters relate to the convolution parameters
%                 corresponding to the convolution function defined by the
%                 'resfun' option (either one or two parameters depending
%                 on function type.
%               - The last parameter is a scale factor for the intensity
%                 If this is omitted, a scale factor of 1 is used;
%
%   kwpars      - A series of 'keywords' and parameters. Specific to this
%                 function is:
%
%               - 'resfun' - determines the convolution / resolution 
%                    function to get S(q,w). It can be either a string: 
%                      'gauss' - gaussian with single fixed (fittable) FWHM
%                      'lor' - lorentzian with single fixed (fittable) FWHM
%                      'voigt' - pseudo-voigt with single fixed (fittable) FWHM
%                      @fun - a function handle satisfying the requirements of
%                             the 'fwhm' parameter of disp2sqw.
%                    NB. For 'gauss' and 'lor' only one fwhm parameter may be
%                    specified. For 'voigt', fwhm = [width lorz_frac]
%                    contains two parameters - the fwhm and lorentzian fraction  
%                    [default: 'gauss']
%               - 'partrans' - a function to transform the fit parameters
%                    This transformation will be applied before each iteration
%                    and the transformed input parameter vector passed to
%                    spinW and the convolution function.
%                    [default: @(y)y  % identity operation]
%
%               In addition, the following parameters are used by this function                         
%                    and will also be passed on to spinw.matparser which will
%                    do the actual modification of spinW model parameters:
%                  
%               - 'mat' - A cell array of labels of spinW named 'matrix' or
%                    matrix elements. E.g. {'J1', 'J2', 'D(3,3)'}. These will
%                    be the model parameters to be varied in a fit, their
%                    order in this cell array will be the same as in the
%                    fit parameters vector.
%                    [default: [] % empty matrix - no model parameters] 
%
%                 All other parameters will be passed to spinW. See the help
%                    for spinw/spinwave, spinw/matparser and spinw/sw_neutron
%                    for more information.
%
%   swobj       The spinwave object which defines the magnetic system to be
%               calculated.
%
% Output:
% -------
%   weight      Array with spectral weight at the q,e points
%               If q and en given:  weight is an nq x ne array, where nq
%                                   is the number of q points, and ne the
%                                   number of energy points
%               If qw given together: weight has the same size and dimensions
%                                     as q{1} i.e. qh
%
% Example:
% --------
%
% tri = sw_model('triAF',[5 1]);                         % J1=5, J2=1 (AFM)
% spinw_pars = {'mat', {'J1', 'J2'}, 'hermit', true, ...
%               'useMex', true, 'optmem', 100};
% [wf,fp] = fit_sqw(w1, @tri.horace_sqw, {[J1 J2 fwhm] spinw_pars});

% Error checking
if ~isnumeric(qh) || ~isnumeric(qk) || ~isnumeric(ql) || ~isnumeric(en)
    error('horace_sqw:BadInput', 'Inputs for qh, qk, ql and en must be numeric vectors.');
elseif numel(qh) ~= numel(qk) || numel(qh) ~= numel(ql) || numel(qh) ~= numel(en)
    error('horace_sqw:BadInput', 'Inputs for qh, qk, ql and en must be same size.');
end

inpForm.fname  = {'resfun' 'partrans' 'mat'};
inpForm.defval = {'gauss'  @(y)y      []};
inpForm.size   = {[1 -2]   [1 1]      [1 -1]};
inpForm.soft   = {false    false      false};

warnState = warning('off','sw_readparam:UnreadInput');
param = sw_readparam(inpForm, varargin{:});

% Sets the number of spinW model parameters. All others are convolution pars.
n_horace_pars = numel(param.mat);
if isempty(n_horace_pars);
    n_horace_pars = 0;
end
model_pars = pars(1:n_horace_pars);
res_pars = pars((n_horace_pars+1):end);
scale_factor = 1;

% Transforms the input parameters
pars = param.partrans(pars);

% Determine which resolution function to use.
if ischar(param.resfun)
    if strncmp(lower(param.resfun), 'gauss', 5) 
        fwhm = res_pars(1);
        if numel(res_pars) > 1
            scale_factor = res_pars(2);
        end
    elseif strncmp(lower(param.resfun), 'lor', 3)
        fwhm = @(emat, cen) lorz_internal(emat, cen, res_pars);
        if numel(res_pars) > 1
            scale_factor = res_pars(2);
        end
    elseif strncmp(lower(param.resfun), 'voigt', 3)
        fwhm = @(emat, cen) voigt_internal(emat, cen, res_pars);
        if numel(res_pars) > 2
            scale_factor = res_pars(3);
        end
    else
        error('horace_sqw:UnknowResFun', ...
            sprintf('Unknown resolution function %s', param.resfun));
    end
else
    fwhm = param.resfun;
    if numel(res_pars) > 0
        scale_factor = res_pars(1);
    end
end

weight = disp2sqw(qh, qk, ql, en, @obj.horace, {model_pars varargin{:}}, fwhm);

if scale_factor ~= 1
    weight = weight * scale_factor;
end

end

%--------------------------------------------------------------------------------------------------
function out = lorz_internal(Emat, center, fwhm)
% Calculates a Lorentzian function for disp2sqw.
    % fwhm should be scalar.
    out = abs(fwhm(1)/pi) ./ (bsxfun(@minus, center, Emat).^2 + fwhm(1)^2);
end


function out = voigt_internal(Emat, center, fwhm)
% Calculates a pseudo-Voigt function for disp2sqw.
    lorfrac = fwhm(2);
    fwhm = fwhm(1);
    sig = fwhm / sqrt(log(256));
    Ediff2 = bsxfun(@minus, center, Emat).^2;
    out = (abs(fwhm/pi) ./ (Ediff2 + fwhm^2)) .* lorfrac + ...
          (exp(-Ediff2 ./ (2*sig^2)) ./ (sig*sqrt(2*pi))) .* (1 - lorfrac);
end
