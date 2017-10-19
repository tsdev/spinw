function y = gauss(x,p)
% normalized Gaussian function
% 
% ### Syntax
% 
% `y = func.gauss(x,p)`
% 
% ### Description
% 
% `y = func.gauss(x,p)` calculates the $y$ values for a Gaussian function
% evaluated at $x$ and with parameters defined in `p`.
% 
% ### Input Arguments
% 
% `x`
% : Coordinate vector where the function will be evaluated.
% 
% `p`
% : Parameter vector with the following elements `p=[I x0 \\sigma]` where:
%   * `I` integrated intensity,
%   * `x0` center,
%   * `\\sigma` standard deviation.
% 
% ### See Also
% 
% [swfunc.pvoigt] \| [swfunc.gaussfwhm]
%

% standard deviation
sigma = p(3);

y = p(1)/sqrt(2*pi)/sigma * exp(-0.5*((x-p(2))/sigma).^2);

end