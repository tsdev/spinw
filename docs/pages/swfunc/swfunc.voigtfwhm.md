---
{title: swfunc.voigtfwhm, link: swfunc.voigtfwhm, summary: normalized function that
    calculates the voigt profile using FWMH values, keywords: sample, sidebar: sw_sidebar,
  permalink: swfunc_voigtfwhm.html, folder: swfunc, mathjax: 'true'}

---

### Syntax

`y = fitfun.voigtfwhm(x,p)`

### Description

The integral of the function is normalized assumind dx = 1.
 

### Input Arguments

`x`
:Input coordinates where the function will be calculated.

`p`
:Parameters:
     A = p(1) integral of the signal assumin dx = 1 (for different
     dx multiply the amplitude with dx to keep the integral
     constant).
     x0 = p(2) peak center positions.
     wG = p(3) FWHM of the Gaussian component.
     wL = p(4) FWHM of the Lorentzian component.

`Conversion`
:sion between different width:

`gamma`
:parameter of the Lorentzian

`gamma`
:= wL/2

`Standard`
:rd deviation of the Gaussian:

`sigma`
:= wG/sqrt(8*ln(2))

### See Also

[func.gauss] \| [func.gaussfwhm]

