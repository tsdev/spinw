function y = gaussfwhm(x,p)
% normalized Gaussian function
% 
% ### Syntax
% 
% `y = func.gaussfwhm(x,p)`
% 
% ### Description
% 
% `y = func.gaussfwhm(x,p)` calculates the $y$ values for a Gaussian
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
% [swfunc.pvoigt] \| [swfunc.gauss]
%

% standard deviation
sigma = p(3)/sqrt(8*log(2));

y = p(1)/sqrt(2*pi)/sigma * exp(-0.5*((x-p(2))/sigma).^2);

end