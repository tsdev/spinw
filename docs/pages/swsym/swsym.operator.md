---
{title: swsym.operator, link: swsym.operator, summary: generates all symmetry elements
    from given space group, keywords: sample, sidebar: sw_sidebar, permalink: swsym_operator,
  folder: swsym, mathjax: 'true'}

---
  
### Syntax
  
`[symOp, symInfo] = swsym.operator(sym)`
  
`[symOp, symInfo] = swsym.operator(sym,fid)`
 
### Description
  
`[symOp, symInfo] = swsym.operator(sym)` generates *all* symmetry
elements from a given set of generators. It also accepts space group
labels or space group index or string of symmetry operators.
  
### Input Arguments
  
`sym`
: Line index in the `symmetry.dat` file or string of the
  symmetry operators or matrix of symmetry generators with dimensions of
  $$[3\times 4\times n_{op}]$$. For example: `sym = 'P n m a'`.
  
`fid`
: If non-zero, the symmetry operators will be printed to the file
  identified by `fid`, the following values are valid:
  * `0`   no printed output (default),
  * `1`   standard output (Command Line),
  * `fid` text file opened before using `fid = fopen(path)`.
  
### Output Arguments
  
`symOp`
: All the symmetry elements in a matrix with dimensions of $$[3\times
  4\times n_{op}]$$.
 
`symInfo`
: Structure that contains additional information about the space 
  group with the following fields:
  * `name`    Name of the space group, if the `swsym.generator`
              function is called with no input, name stores the name of
              all space groups from `symmetry.dat` file in a cell.
  * `str`     The string of the symmetry operations.
  * `num`     The line index in the `symmetry.dat` file.
  
### See Also
  
[swsym.generator](swsym_generator)
 

{% include links.html %}
