---
{title: swsymposition( ), keywords: sample, summary: generates symmetry equivalent
    positions, sidebar: sw_sidebar, permalink: +swsym_position.html, folder: +swsym,
  mathjax: 'true'}

---
generates symmetry equivalent positions
 
[r, aIdx, opInfo] = SWSYM.POSITION(sym, r0, fid, tol)
  
It generates all symmetry equivalent atomic positions from a given
symmetry number and coordinates of the input atoms. If fid is defined,
the result is printed onto the command window.
 
Input:
 
symOp         Matrix containing the symmetry operators 
r0            Atomic position in lattice units, dimensions are [3 nAtom].
fid           Optional input, the file identifier to print the result.
                  0   No output printed (Default)
                  1   Output printed onto the screen (Command Window)
                  fid Use with the following command: fid = fopen(path)
tol           Tolerance, distance within two atoms are considered
              identical, default is 1e-5 lattice unit. Necessary for
              badly defined atomic positions (when atoms are not exactly
              on the symmetry element) and to avoid numerical errors.
 
Output:
 
rSym          All generated atomic positions, dimensions are
              [3 nGenAtom].
aIdx          The index of the symmetry inequivalent position for every
              generated position, dimensions are [1 nGenAtom].
opInfo        Structure with the following fields:
  ismoved         Cell, where each element contains a vector with logical
                  value, whether the given operator moved the atom or
                  not. Each vector has a dimensions of [1 nSym], where
                  the nSym is multiplicity of the general position.
  opmove          The rotation operators for every moved atom, dimensions
                  are [3 3 nGenAtom].
 
See also SPINW, SWSYM.OPERATOR, SPINW.ATOM, SPINW.MATOM.
 
