---
{title: swfunc.voigtfwhm( ), summary: normalized function that calculates the voigt
    profile using FWMH values, keywords: sample, sidebar: sw_sidebar, permalink: swfunc_voigtfwhm.html,
  folder: swfunc, mathjax: 'true'}

---
normalized function that calculates the voigt profile using FWMH values
 
y = FITFUN.VOIGTFWHM(x,p)
 
The integral of the function is normalized assumind dx = 1.
 
Input:
 
x     Input coordinates where the function will be calculated.
p     Parameters:
          A = p(1) integral of the signal assumin dx = 1 (for different
          dx multiply the amplitude with dx to keep the integral
          constant).
          x0 = p(2) peak center positions.
          wG = p(3) FWHM of the Gaussian component.
          wL = p(4) FWHM of the Lorentzian component.
 
Conversion between different width:
gamma parameter of the Lorentzian
gamma = wL/2
Standard deviation of the Gaussian:
sigma = wG/sqrt(8*ln(2))
 
See also FUNC.GAUSS, FUNC.GAUSSFWHM.

