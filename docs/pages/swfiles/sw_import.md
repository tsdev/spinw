---
{title: sw_import, link: sw_import, summary: create SpinW object from .cif and FullProf
    Studio .fst files, keywords: sample, sidebar: sw_sidebar, permalink: sw_import,
  folder: swfiles, mathjax: 'true'}

---
  
### Syntax
  
`obj = sw_import(fname, {toplot})`
  
### Description
  
`obj = sw_import(fname, {toplot})` can import Crystallographic
Information Framework (.cif) files or FullProf Studio (.fst) files. It is
also able to read .cif files from a web link. 
  
### Input Arguments
  
`fName`
: String that contains the location of the source file, e.g. the full
  path to the file or a web address.
  
`toPlot`
: If `true` the imported structure will be plotted, default value is
  `false`.
 

{% include links.html %}
