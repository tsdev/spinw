---
{title: swsym.generator, link: swsym.generator, summary: returns symmetry operators
    of a given space group, keywords: sample, sidebar: sw_sidebar, permalink: swsym_generator,
  folder: swsym, mathjax: 'true'}

---
  
### Syntax
  
`[symOp, symInfo] = swsym.generator(sym)`
  
`[symOp, symInfo] = swsym.generator(sym,fid)`
 
### Description
  
`[symOp, symInfo] = swsym.generator(sym)` gives the symmetry operators
based on a given space group number or a string of symmetry operators.
Without arguments, the function returns the name of all space groups
stored in `symmetry.dat` file.
   
`[symOp, symInfo] = swsym.generator(sym,fid)` also prints the symmetry
operators to the file identified by `fid`.
 
### Input Arguments
  
`sym`
: Either the label of the space group or the index from
  the [International Tables of Crystallography](http://it.iucr.org/A/) or
  string containing the space group operators in the same format as used
  in the `symmetry.dat` file (for details see [swsym.str](swsym_str)).
  
`fid`
: If non-zero, the symmetry operators will be printed to the file
  identified by `fid`, the following values are valid:
  * `0`   no printed output (default),
  * `1`   standard output (Command Line),
  * `fid` text file opened before using `fid = fopen(path)`.
  
### Output Arguments
  
`symOp`
: Symmetry operators in a matrix with dimensions of $$[3\times 4\times
  n_{op}]$$.
 
`symInfo`
: Structure that contains additional information about the space 
  group with the following fields:
  * `name`    Name of the space group, if the `swsym.generator`
              function is called with no input, name stores the name of
              all space groups from `symmetry.dat` file in a cell.
  * `str`     The string of the symmetry operations.
  * `num`     The line index in the `symmetry.dat` file.
  
### See Also
  
[swsym.add](swsym_add) \| [spinw](spinw) \| [spinw.gencoupling](spinw_gencoupling) \| [swsym.position](swsym_position)
 

{% include links.html %}
