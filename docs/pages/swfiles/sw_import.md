---
{title: sw_import( ), summary: create SpinW object from .cif and FullProf Studio .fst
    files, keywords: sample, sidebar: sw_sidebar, permalink: swfiles_sw_import.html,
  folder: swfiles, mathjax: 'true'}

---
create SpinW object from .cif and FullProf Studio .fst files
 
At present the function can import Crystallographic Information Framework
(.cif) files or FullProf Studio (.fst) files. It is able to read .cif
files from a web link.
 
obj = SW_IMPORT(fName, {toPlot})
 
Input:
 
fName     String, contain the file location and name of the .fst file.
toPlot    If true the structure will be plotted, default is false.
 
