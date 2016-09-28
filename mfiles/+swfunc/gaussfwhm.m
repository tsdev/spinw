function y = gaussfwhm(x,p)
% normalized gaussian function defined by the FWHM
%
% y = FUNC.GAUSSFWHM(x,p)
%
% The integral of the function is normalized assumind dx = 1.
%
% Input:
% x         Coordinate vector where the function will be evaluated.
% p         Parameter vector: [I Centre FWHM].
%
% See also SWFUNC.VOIGT, SWFUNC.GAUSS.

% standard deviation
sigma = p(3)/sqrt(8*log(2));

y = p(1)/sqrt(2*pi)/sigma * exp(-0.5*((x-p(2))/sigma).^2);

end