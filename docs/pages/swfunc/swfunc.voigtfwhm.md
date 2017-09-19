---
{title: swfunc.voigtfwhm, link: swfunc.voigtfwhm, summary: normalized function that
    calculates the voigt profile using FWMH values, keywords: sample, sidebar: sw_sidebar,
  permalink: swfunc_voigtfwhm.html, folder: swfunc, mathjax: 'true'}

---

### Syntax

` `

### Description

 

### Input Arguments

% `x`
:  Input coordinates where the function will be calculated.

% `p`
:  Parameters:

% `A`
:p(1) integral of the signal assumin dx = 1 (for different

% `dx`
:ultiply the amplitude with dx to keep the integral

% ``
:tant).

% `x0`
: p(2) peak center positions.

% `wG`
: p(3) FWHM of the Gaussian component.

% `wL`
: p(4) FWHM of the Lorentzian component.

% `Conversion`
:ersion between different width:

% `gamma`
:a parameter of the Lorentzian

% `gamma`
:a = wL/2

% `Standard`
:dard deviation of the Gaussian:

% `sigma`
:a = wG/sqrt(8*ln(2))

### See Also

[func.gauss] and [func.gaussfwhm]

