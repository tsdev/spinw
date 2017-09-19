---
{title: swsym.operator, link: swsym.operator, summary: calculates all symmetry operators
    or general positions for a space group, keywords: sample, sidebar: sw_sidebar,
  permalink: swsym_operator.html, folder: swsym, mathjax: 'true'}

---

### Syntax

` `

### Description



### Input Arguments

% `sym`
: Line index in the symmetry.dat file or string of the
 ymmetry operators or matrix of symmetry generators with
 imensions of [3 4 nOp].
 or example:
    sym = 'P b n m';

% `fid`
: For printing the symmetry operators:
    0   no printed output (Default)
    1   standard output (Command Line)
    fid text file opened before with the fid = fopen(path)

### Output Arguments

symOp         The rotational part of the symmetry operators, dimensions
            are [3 3 nSym].
symInfo       Structure containing additional information about the space
            group with the following fields:
name            Name of the space group in string. If function called
                with no input, name stores the name of all spase groups
                from symmetry.dat in a cell.
str             The string of the symmetry operations.
num             The index of the symmetry in the symmetry.dat file.

### See Also

[spinw](spinw.html), [swsym.position](swsym_position.html), [spinw.atom](spinw_atom.html), [spinw.matom](spinw_matom.html) and [swsym.generator](swsym_generator.html)
SWSYM.POINT.

