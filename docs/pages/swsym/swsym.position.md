---
{title: swsym.position, link: swsym.position, summary: generates symmetry equivalent
    positions, keywords: sample, sidebar: sw_sidebar, permalink: swsym_position, folder: swsym,
  mathjax: 'true'}

---
  
### Syntax
  
`[r, aIdx, opInfo] = swsym.position(sym,r0)`
  
`[r, aIdx, opInfo] = swsym.position(sym,r0,fid)`
 
`[r, aIdx, opInfo] = swsym.position(sym,r0,fid,tol)`
 
### Description
  
`[r, aIdx, opInfo] = swsym.position(sym, r0, fid, tol)` generates all
symmetry equivalent atomic positions from a given space group and
coordinates of the symmetry inequivalent atoms. If `fid` is defined, the
result are printed onto the corresponding file.
  
### Input Arguments
  
`sym`
: Either the label of the space group or the index from
  the [International Tables of Crystallography](http://it.iucr.org/A/) or
  string containing the space group operators in the same format as used
  in the `symmetry.dat` file (for details see [swsym.str](swsym_str)).
  
`r0`
: Atomic position in lattice units in a matrix with dimensions of
  $$[3\times n_{atom}]$$.
  
`fid`
: If non-zero, the symmetry operators will be printed to the file
  identified by `fid`, the following values are valid:
  * `0`   no printed output (default),
  * `1`   standard output (Command Line),
  * `fid` text file opened before using `fid = fopen(path)`.
  
`tol`
: Tolerance, distance within two atoms are considered
  identical, default value is $$10^{-5}$$ lattice unit. Necessary to check
  for badly defined atomic positions (when atoms are not exactly on the
  symmetry element) and to avoid numerical errors.
  
### Output Arguments
  
`r`
: All generated atomic positions stored in a matrix with dimensions of
  $$[3\times n_{genAtom}]$$.
 
`aIdx`
: The index of the symmetry inequivalent position for every
  generated position, stored in a row vector with $$n_{genAtom}$$ number of
  elements.
 
`opInfo`
: Structure with the following fields:
  * `ismoved`     Cell, where each element contains a vector with logical
                  values, whether the given operator moved the atom or
                  not. Each vector has a dimensions of $$[1\times n_{sym}]$$, where
                  the $$n_{sym}$$ is multiplicity of the general position.
  * `opmove`      The rotation operator that moved the original atom the
                  equivalent position stored in a matrix with dimensions
                  of $$[3\times 3\times n_{genAtom}]$$.
  
### See Also
  
[swsym.operator](swsym_operator)
 

{% include links.html %}
