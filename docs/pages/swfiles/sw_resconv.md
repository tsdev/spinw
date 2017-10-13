---
{title: sw_resconv, link: sw_resconv, summary: convolution of a matrix and a Gaussian,
  keywords: sample, sidebar: sw_sidebar, permalink: sw_resconv, folder: swfiles, mathjax: 'true'}

---
  
### Syntax
  
`Mout = sw_resconv(M,x,dx,func)`
  
### Description
  
`Mout = sw_resconv(M,x,dx,func)` convolutes a 2D matrix with a Gaussian
along the first dimension of the matrix. The convolution keeps the
integrated intensity $$\sum I\cdot dx$$ constant. It assumes the `x` vector
contains the center points of the bins and the distances between the
generated bin edges is calculated by interpolating from the distances
between the given `x` bin center positions.
  
### Input Arguments
  
`M`
: Arbitrary matrix with dimensions of $$[m_1\times m_2]$$.
  
`x`
: Column vector of coordinates along the first dimension of the
  matrix.
  
`dx`
: FWHM value of the Gaussian as a
  function of $$dx$$. Either a function handle with a header `fwhm =
  dx(xVal)` or a vector of polynomial coefficients that produces the
  right standard deviation. In this case in the function the following
  line will be executed `fwhm = polyval(dx,xVal)` or a constant FWHM
  value.
 
  The standard deviation of the Gaussian is calculated from the given
  $$dx$$ value using the formula $$\sigma_G = fwhm_G/2/\sqrt{2\cdot log(2)}
  \sim fwhm_G\cdot 0.424661$$
  If a general resolution function is provided in the `func` argument,
  it will be called as `y = func(x,[1 x0 fwhm])`. In this case the `fwhm`
  can be a row vector and the meaning of the different parameters will
  depend on `func`.
  
`func`
: Resolution function shape with header `y = func(x,p)`,
  where `x` is a column vector, `p` is a row vector of parameters. The
  meaning of the first 2 elements of the parameter vector are fixed.
 
  * `p(1)` integral of the function.
  * `p(2)` center of mass position of the function.
 
  Optional, default value is `@swfunc.gaussfwhm`.
  
### Output Arguments
  
`Mout`
: Matrix with same dimensions as the input `M`.
  
### See Also
  
[sw_res](sw_res) \| [swfunc.gaussfwhm](swfunc_gaussfwhm)
 
*[FWHM]: Full Width at Half Maximum
 

{% include links.html %}
