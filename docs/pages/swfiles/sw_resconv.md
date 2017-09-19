---
{title: sw_resconv( ), link: sw_resconv, summary: Convolute Gaussian with variable
    width along the first dimension of a matrix, keywords: sample, sidebar: sw_sidebar,
  permalink: sw_resconv.html, folder: swfiles, mathjax: 'true'}

---

### Syntax

` `

### Description

assumes the x vector contains the center points of the bins and the
distances between the generated bin edges is calculated by interpolating
from the distances between the given x bin center positions.
 

### Input Arguments

% `M`
: Arbitrary matrix with dimensions of (m1,m2).

% `x`
: Column vector of coordinates along the first dimension of the

% ``
:x.

% `dx`
: Full width at half maximum (FWHM) value of the Gaussian as a

% `function`
:ion of dx. Either a function handle with a header:
 whm = dx(xVal)

% `or`
:vector of polynomial coefficients that produces the right

% `standard`
:ard deviation. In this case in the function the following line

% `will`
:be executed:
 whm = polyval(dx,xVal)

% `or`
:constant fwhm value.

% `The`
:tandard deviation of the Gaussian is calculated from the given

% `dx`
:lue using the following formula:
 tdG = fwhmG/2/sqrt(2*log(2)) ~ fwhmG/2.35482

% `If`
:general resolution function is provided in the func argument,

% `it`
:ll be called as:
  = func(x,[1 x0 fwhm]);

% `In`
:is case the fwhm can be a row vector and the meaning of the

% `different`
:rent parameters will depend on func.

% `func`
: Resolution function shape with header:
  = func(x,p)

% `Where`
: x is a column vector, p is a row vector of parameters. The

% `meaning`
:ng of the first 2 elements of the parameter vector are fixed.
 (1)    integral of the function.
 (2)    center of mass position of the function.

% `Optional,`
:nal, default value is @swfunc.gaussfwhm.

### Output Arguments

M     Matrix with same dimensions as the input storing the convoluted
data.

### See Also

[sw_res](sw_res.html)

