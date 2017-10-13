---
{title: swfunc.gaussfwhm, link: swfunc.gaussfwhm, summary: normalized Gaussian function,
  keywords: sample, sidebar: sw_sidebar, permalink: swfunc_gaussfwhm, folder: swfunc,
  mathjax: 'true'}

---
  
### Syntax
  
`y = func.gaussfwhm(x,p)`
  
### Description
  
`y = func.gaussfwhm(x,p)` calculates the $$y$$ values for a Gaussian
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
  
[swfunc.pvoigt](swfunc_pvoigt) \| [swfunc.gauss](swfunc_gauss)
 

{% include links.html %}
