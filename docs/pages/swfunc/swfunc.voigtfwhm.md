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
:nput coordinates where the function will be calculated.

`p`
:arameters:
    A = p(1) integral of the signal assumin dx = 1 (for different
    dx multiply the amplitude with dx to keep the integral
    constant).
    x0 = p(2) peak center positions.
    wG = p(3) FWHM of the Gaussian component.
    wL = p(4) FWHM of the Lorentzian component.

`Conversion`
:ion between different width:

`γ`
: parameter of the Lorentzian

`γ`
: = wL/2

`Standard`
:d deviation of the Gaussian:

`sigma`
: wG/sqrt(8*ln(2))

### See Also

[func.gauss] \| [func.gaussfwhm]

