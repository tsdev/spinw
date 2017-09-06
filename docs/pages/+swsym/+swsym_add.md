---
{title: swsymadd( ), keywords: sample, summary: saves user defined symmetry operators,
  sidebar: sw_sidebar, permalink: +swsym_add.html, folder: +swsym, mathjax: 'true'}

---
saves user defined symmetry operators
 
sym = SWSYM.ADD(symStr, {symName})
 
It saves the symmetry generators in symStr into the symmetry.dat file and
returns the line number of the space group in the symmetry.dat file.
 
Input:
symStr        String, that contains the operators of the space group. If
              not only the generators are given, a possible set of
              generators will be determined and only those will be saved.
symName       Label for the space group.
 
See also SWSYM.GENERATOR, SWSYM.GENREDUCE.
 
