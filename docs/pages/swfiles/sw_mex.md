---
{title: sw_mex, link: sw_mex, summary: compiles and tests the mex files, keywords: sample,
  sidebar: sw_sidebar, permalink: sw_mex, folder: swfiles, mathjax: 'true'}

---
  
### Syntax
  
`sw_mex(Name,Value)`
  
### Description
  
`sw_mex(Name,Value)` compiles and tests the generated mex files. The
compiled mex files will speed up the [spinw.spinwave](spinw_spinwave) function. The
expected speedup is larger for smaller magnetic unit cells. Once the mex
files are compiled, use the `swpref.setpref('usemex',true)` command to
switch on using mex files in [spinw.spinwave](spinw_spinwave).
  
### Name-Value Pair Arguments
  
`'test'`
: If `true`, the compiled .mex files will be tested. Default is
  `false`.
  
`'swtest'`
: If `true`, 3 spin wave calculation will run with and without .mex
  files and the results will be compared. Default is `false`.
  
### See Also
  
[swpref]
 

{% include links.html %}
