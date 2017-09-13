---
{title: sw_intsf( ), summary: integrates the structure factor along given Q directions,
  keywords: sample, sidebar: sw_sidebar, permalink: sw_intsf.html, folder: swfiles,
  mathjax: 'true'}

---
integrates the structure factor along given Q directions
 
sFact = SW_INTSF(sFact, 'Option1', Value1, ...) 
 
Options:
 
axis      The index of the axis, along which the data is summed within
          the given range. {1, 2, 3} is for {a*, b*, c*} respectively.
          Zero is to powder average the data. Default is 1. More than one
          axis can be given.
range     Data range in r.l.u. along selected dimension to integrate,
          default is the full data range.
type      Type of data to be used for the integration.
              'sf'    integrate the calculated structure factor. (default)
              'perp'  integrate the perpendicular component of the structure
                      factor to Q.
 
Output:
          sFact contains the following new fields:
          int     The integrated data in a matrix.
          hklint  A cell containing the three axis vectors in r.l.u. and
 
Example:
  sFact = cryst.structfact;
  sFact = sw_intsf(sFact,'axis',[1 2],'range',[0 1; 0 2]);
 
  The above code will integrate the structure factor along a* between 0
  and 1, along b* between 0 and 2.
 
See also SPINW, SPINW.STRUCTFACT, SW_PLOTSF.
 
