function y = lorfwhm(x,p)
% normalized Lorentzian function
% 
% ### Syntax
% 
% `y = func.lorfwhm(x,p)`
% 
% ### Description
% 
% `y = func.lorfwhm(x,p)` calculates the $y$ values for a Lorentzian
% function evaluated at $x$ and with parameters defined in `p`.
% 
% ### Input Arguments
% 
% `x`
% : Coordinate vector where the function will be evaluated.
% 
% `p`
% : Parameter vector with the following elements `p=[I x0 FWHM]` where:
%   * `I`       integrated intensity,
%   * `x0`      center,
%   * `FWHM`    Full Width at Half Maximum value.
% 
% ### See Also
% 
% [swfunc.pvoigt] \| [swfunc.gaussfwhm]
%

y = p(1)/(pi*p(3))./(1+((x-p(2))/p(3)).^2);

end