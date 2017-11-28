function y = pvoigt(x,p)
% pseudovoigt function
% 
% ### Syntax
% 
% `y = func.pvoigt(x,p)`
% 
% ### Description
% 
% `y = func.pvoigt(x,p)` calculates the $y$ values for a pseudovoigt
% function evaluated at $x$ and with parameters defined in `p`. The
% Gaussian and Lorentzian functions are normalized to amplitude 1 before
% the mixing.
% 
% ### Input Arguments
% 
% `x`
% : Coordinate vector where the function will be evaluated.
% 
% `p`
% : Parameter vector with the following elements `p=[A x0 wG wL mu]` where:
%   * `A`       amplitude of the signal,
%   * `x0`      center of the peak,
%   * `wG`      FWHM value of the Gaussian, has to be positive,
%   * `wL`      FWHM value of the Lorentzian, has to be positive,
%   * `mu`      mixing constant, `mu=1` for pure Lorenzian, `mu=0` for pure
%               Gaussian, the value has to be within the $(0,1)$ range.
% 
% ### See Also
% 
% [swfunc.gaussfwhm] \| [swfunc.lorfwhm] \| [swfunc.voigtfwhm]
%
% *[FWHM]: Full Width at Half Maximum
%

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