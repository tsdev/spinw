---
{title: swfunc.lorfwhm, link: swfunc.lorfwhm, summary: normalized Lorentzian function,
  keywords: sample, sidebar: sw_sidebar, permalink: swfunc_lorfwhm, folder: swfunc,
  mathjax: 'true'}

---
  
### Syntax
  
`y = func.lorfwhm(x,p)`
  
### Description
  
`y = func.lorfwhm(x,p)` calculates the $$y$$ values for a Lorentzian
function evaluated at $$x$$ and with parameters defined in `p`.
  
### Input Arguments
  
`x`
: Coordinate vector where the function will be evaluated.
  
`p`
: Parameter vector with the following elements `p=[I x0 FWHM]` where:
  * `I`       integrated intensity,
  * `x0`      center,
  * `FWHM`    Full Width at Half Maximum value.
  
### See Also
  
[swfunc.pvoigt](swfunc_pvoigt) \| [swfunc.gaussfwhm](swfunc_gaussfwhm)
 

{% include links.html %}
