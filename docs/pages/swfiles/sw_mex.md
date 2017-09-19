---
{title: sw_mex( ), link: sw_mex, summary: compiles the mex files and test them, keywords: sample,
  sidebar: sw_sidebar, permalink: sw_mex.html, folder: swfiles, mathjax: 'true'}

---

### Syntax

` `

### Description

expected speedup is larger for smaller magnetic unit cells. Once the mex
files are compiled, use the swpref.setpref('usemex',true) command to
switch to mex in spinw.spinwave.
 

### Name-Value Pair Arguments

% `test`
:  If true, the compiled .mex files will be tested. Default is

% ``
:

% `swtest`
:  If true, 3 spin wave calculation will run with and without .mex

% `files`
:d the results will be compared. Default is false.

### See Also

[swpref]

