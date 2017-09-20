---
{title: sw_filelist( ), link: sw_filelist, summary: lists spinw data in the Matlab
    workspace or in a .mat file, keywords: sample, sidebar: sw_sidebar, permalink: sw_filelist.html,
  folder: swfiles, mathjax: 'true'}

---

### Syntax

`list = sw_filelist({Name,Value)`

### Description

 

### Name-Value Pair Arguments

`fName`
: To check data stored in a .mat file, fName contains the path as
  a string.

`sort`
: To sort according to which column (positive ascending, negative
  descending):
      +/-1    variable name,
      +/-2    title,
      +/-3    creation date,
      +/-4    completion date, default.

### Output Arguments

list      Cell of strings, lists each simulation data in the Matlab
          memory, even data stored in cells.

### See Also

[spinw.anneal](spinw_anneal.html) \| [spinw.spinwave](spinw_spinwave.html)

