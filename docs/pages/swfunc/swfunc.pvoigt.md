---
{title: swfunc.pvoigt, link: swfunc.pvoigt, summary: pseudovoigt function, keywords: sample,
  sidebar: sw_sidebar, permalink: swfunc_pvoigt, folder: swfunc, mathjax: 'true'}

---
  
### Syntax
  
`y = func.pvoigt(x,p)`
  
### Description
  
`y = func.pvoigt(x,p)` calculates the $$y$$ values for a pseudovoigt
function evaluated at $$x$$ and with parameters defined in `p`. The
Gaussian and Lorentzian functions are normalized to amplitude 1 before
the mixing.
  
### Input Arguments
  
`x`
: Coordinate vector where the function will be evaluated.
  
`p`
: Parameter vector with the following elements `p=[A x0 wG wL mu]` where:
  * `A`       amplitude of the signal,
  * `x0`      center of the peak,
  * `wG`      FWHM value of the Gaussian, has to be positive,
  * `wL`      FWHM value of the Lorentzian, has to be positive,
  * `mu`      mixing constant, `mu=1` for pure Lorenzian, `mu=0` for pure
              Gaussian, the value has to be within the $$(0,1)$$ range.
  
### See Also
  
[swfunc.gaussfwhm](swfunc_gaussfwhm) \| [swfunc.lorfwhm](swfunc_lorfwhm) \| [swfunc.voigtfwhm](swfunc_voigtfwhm)
 
*[FWHM]: Full Width at Half Maximum
 

{% include links.html %}
