function y = lorfwhm(x,p)
% normalized Lorentzian using FWHM
%
% y = LORFWHM(x,p)
%
% The function is normalized to integral 1 for dx=1.
%
% Input:
%
% x     Vector of coordinate values where the function is evaluated.
% p     Parameter vector with values  p = [I Centre FWHM].
%
%

y = p(1)/(pi*p(3))./(1+((x-p(2))/p(3)).^2);

end