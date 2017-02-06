function y = pvoigt(x,p)
% pseudovoigt function
%
% y = pvoigt(x,p)
%
% The Gaussian and Lorentzian functions are normalized to amplitude 1
% before the mixing.
%
% Input parameters:
%
% x     Coordinate values.
% p 	Function parameter values: p = [A x0 wG wL mu], where:
%           A       Amplitude of the signal.
%           x0      Center of the peak.
%           wG      FWHM value of the Gaussian, has to be positive.
%           wL      Width of the Lorentzian has to be positive.
%           mu      Mixing constant, mu = 1 for pure Lorentzian, mu = 0 for
%                   pure Gaussian. It has to be within the range of [0 1].
%
% See also SWFUNC.GAUSS, SWFUNC.LORFWHM.

A  = p(1);
x0 = p(2);
wG = p(3);
wL = p(4);
mu = p(5);

if mu < 0 || mu > 1
    error('Mixing constant has to be within the [0 1] range!')
end

if wG < 0 || wL < 0
    error('The width of Lorentzian and Gaussian has to be positive!')
end

yV = func.gaussfwhm(x,[1 x0 wG]);
yL = func.lorfwhm(x,[1 x0 wL]);
y  =  A*(mu*yL + (1-mu)*yV);

end