---
{title: swsym.generator( ), summary: returns symmetry operators of a given space group,
  keywords: sample, sidebar: sw_sidebar, permalink: swsym_generator.html, folder: +swsym,
  mathjax: 'true'}

---
returns symmetry operators of a given space group
 
[symOp, symInfo] = SWSYM.GENERATOR(sym, {fid})
 
It gives the symmetry elements based on the space group number or given
list of symmetry operators. Without arguments, returns the name of all
space groups stored in symmetry.dat file.
 
Input:
 
sym           Either the label of the space group or the index from
              the International Tables of Crystallography or string
              containing the space group operators in the same format as
              used in the symmetry.dat file.
fid           For printing the symmetry operators:
                  0   no printed output (Default)
                  1   standard output (Command Line)
                  fid text file opened before with the fid = fopen(path)
 
Output:
 
symOp         Matrices defining the symmetry operators, dimensions are 
              [3 4 nOp].
symInfo       Structure containing additional information about the space 
              group with the following fields:
  name            Name of the space group in string. If function called
                  with no input, name stores the name of all spase groups
                  from symmetry.dat in a cell.
  str             The string of the symmetry operations.
  num             The index of the symmetry in the symmetry.dat file.
 
See also SWSYM.ADD, SPINW, SPINW.GENCOUPLING, SWSYM.POSITION.
 
