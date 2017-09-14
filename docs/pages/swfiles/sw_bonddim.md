---
{title: sw_bonddim( ), summary: find dimensionality of a periodic bond network, keywords: sample,
  sidebar: sw_sidebar, permalink: sw_bonddim.html, folder: swfiles, mathjax: 'true'}

---
 
L = SW_BONDDIM(C, {nAtom})
 
The function splits the given periodic bond network into disjunct
subsystems and determines the dimensionality of each subsystem.
 
Input:
 
C         Bond list in a matrix with dimensions [5 nBond]. The meaning of
          the rows:
              #1:#3   Lattice translations between the coupled atoms in
                      lattice units (always integer).
              #4      Index of the bond starting atom.
              #5      Index of the bond end atom.
          For example for a chain along b-axis on a Bravais lattice:
              C = [1;1;0;1;0]
nAtom     Number of atoms in the unit cell. If not given, the maximum
          atom index from the bond list is taken.
 
Output:
 
L         Struct with the number of elements equal to the number of
          subsystems, it has the following fields:
              D       Dimensionality of the subsystem (0<=D<=3).
              base    Basis vectors spanning the subsystem stored in a
                      [3 D] matrix where each column denotes a basis
                      vector.
              site    List of sites that belong to the subsystem.
 

