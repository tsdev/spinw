---
{title: swfunc.pvoigt( ), summary: pseudovoigt function, keywords: sample, sidebar: sw_sidebar,
  permalink: swfunc_pvoigt.html, folder: swfunc, mathjax: 'true'}

---
pseudovoigt function
 
y = pvoigt(x,p)
 
The Gaussian and Lorentzian functions are normalized to amplitude 1
before the mixing.
 
Input parameters:
 
x     Coordinate values.
p 	Function parameter values: p = [A x0 wG wL mu], where:
          A       Amplitude of the signal.
          x0      Center of the peak.
          wG      FWHM value of the Gaussian, has to be positive.
          wL      Width of the Lorentzian has to be positive.
          mu      Mixing constant, mu = 1 for pure Lorentzian, mu = 0 for
                  pure Gaussian. It has to be within the range of [0 1].
 
See also SWFUNC.GAUSS, SWFUNC.LORFWHM.
