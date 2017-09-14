---
{title: sw_idata( ), summary: creates iData object, keywords: sample, sidebar: sw_sidebar,
  permalink: sw_idata.html, folder: swfiles, mathjax: 'true'}

---
 
[omega, swConv] = SW_IDATA(spectrum, 'option1', value1 ...) 
 
It creates iData object (<a href=http://ifit.mccode.org>http://ifit.mccode.org</a>) 
and convolutes the spectra with a fixed instrumental resolution, assuming
the energy and Q axis are linear.
 
Input:
 
spectrum      Calculated spin wave spectrum, struct type.
 
Options:
 
fwhmE         Full width half maximum of the Gaussian energy
              resolution. Works properly only for linear energy axis.
fwhmQ         Full width half maximum of the Gaussian momentum
              transfer resolution in A^-1 units. Works properly only
              for linear scans in reciprocal space. Be carefull for
              example for scans like [1, QK, 0] where the equal QK
              steps give unequal steps in the A^-1 reciprocal space.
nInterp       Number of axis subdivision before convolution, equal
              for Q and E.
 

