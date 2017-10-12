function y = gauss(x,p)
% normalized gaussian function
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
%
%   * `I` Integrated intensity.
%   * `x0` Center.
%   * `\\sigma` Standard deviation.
% 
% ### See Also
% 
% [swfunc.voigt] \| [swfunc.gaussfwhm]
%

% normalized gaussian function
%
% y = FUNC.GAUSS(x,p)
%
% The integral of the function is normalized assuming dx = 1.
%
% Input:
% x         Coordinate vector where the function will be evaluated.
% p         Parameter vector: [I Centre sigma].
%
% See also SWFUNC.VOIGT, SWFUNC.GAUSSFWHM.

% standard deviation
sigma = p(3);

y = p(1)/sqrt(2*pi)/sigma * exp(-0.5*((x-p(2))/sigma).^2);

end