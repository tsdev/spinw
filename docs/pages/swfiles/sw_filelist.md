---
{title: sw_filelist( ), summary: lists spinw data in the Matlab workspace or in a
    .mat file, keywords: sample, sidebar: sw_sidebar, permalink: sw_filelist.html,
  folder: swfiles, mathjax: 'true'}

---
lists spinw data in the Matlab workspace or in a .mat file
 
list = SW_FILELIST({'Option1',Value1,...)
 
 
Options:
 
fName     To check data stored in a .mat file, fName contains the path as
          a string.
sort      To sort according to which column (positive ascending, negative
          descending):
              +/-1    variable name,
              +/-2    title,
              +/-3    creation date,
              +/-4    completion date, default.
 
 
Output:
 
list      Cell of strings, lists each simulation data in the Matlab
          memory, even data stored in cells.
 
See also SPINW.ANNEAL, SPINW.SPINWAVE.
 
