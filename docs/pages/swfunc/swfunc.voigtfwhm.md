---
{title: swfunc.voigtfwhm, link: swfunc.voigtfwhm, summary: normalized voigt function,
  keywords: sample, sidebar: sw_sidebar, permalink: swfunc_voigtfwhm.html, folder: swfunc,
  mathjax: 'true'}

---
 
### Syntax
 
`y = fitfun.voigtfwhm(x,p)`
 
### Description
 
`y = fitfun.voigtfwhm(x,p)` calculates the voigt function. The width
parameters define the FWHM value. The
integral of the function is normalized assuming that $$dx = 1$$. The
conversion between different width:
 
* gamma parameter of the Lorentzian $$\gamma = w_L/2$$
* standard deviation of the Gaussian $$\sigma = w_G/\sqrt{8\cdot\ln(2)}$$
 
### Input Arguments
 
`x`
: Input coordinates where the function will be calculated.
 
`p`
: Parameters in a vector with elements `[A x0 wG wL]`:
 
  * `A` integral of the output assuming $$dx=1$$ (for different $$dx$$
     multiply the amplitude with $$dx$$ to keep the integral constant).
  * `x0` peak center positions.
  * `wG` FWHM of the Gaussian component.
  * `wL` FWHM of the Lorentzian component.
 
 
### See also
 
[swfunc.gauss](swfunc_gauss.html) \| [swfunc.gaussfwhm](swfunc_gaussfwhm.html)
 
*[FWHM]: Full Width at Half Maximum
 

