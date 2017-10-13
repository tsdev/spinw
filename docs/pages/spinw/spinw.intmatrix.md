---
{title: spinw.intmatrix method, link: spinw.intmatrix, summary: generates interaction
    matrix, keywords: sample, sidebar: sw_sidebar, permalink: spinw_intmatrix, folder: spinw,
  mathjax: 'true'}

---
  
### Syntax
  
`[SS, SI, RR] = intmatrix(obj,Name,Value)`
  
### Description
  
`[SS, SI, RR] = intmatrix(obj,Name,Value)` lists the bonds and generates
the corresponding exchange matrices by applying the bond symmetry
operators on the stored matrices. Also applies symmetry on the single ion
anisotropies and can generate the representation of bonds, anistropies
and atomic positions in an arbitrary supercell. The output argument `SS`
contains the different types of exchange interactions separated into
different fields, such as `SS.DM` for the Dzyaloshinskii-Moriya
interaction, `SS.iso` for Heisenberg exchange, `SS.all` for general
exchange etc.
  
### Input Arguments
  
`obj`
: [spinw](spinw) object.
  
### Name-Value Pair Arguments
  
`'fitmode'`
: Can be used to speed up calculation, modes:
  * `true`    No speedup (default).
  * `false`   For the interactions stored in `SS`, only the
              `SS.all` field is calculated.
  
`'plotmode'`
: If `true`, additional rows are added to `SS.all`, to identify
  the couplings for plotting. Default is `false`.
  
`'sortDM'`
: If true each coupling is sorted for consistent plotting of
  the DM interaction. Sorting is based on the `dR` bond vector that
  points from `atom1` to `atom2`, for details see [spinw.coupling](spinw_coupling).
  After sorting `dR` vector components fulfill the following rules in
  hierarchical order:
  1. `dR(x) > 0`
  2. `dR(y) > 0`
  3. `dR(z) > 0`.
 
  Default is `false`.
  
`'zeroC'`
: Whether to output bonds with assigned matrices that have only
  zero values. Default is `false`.
  
`'extend'`
: If `true`, all bonds in the magnetic supercell will be
  generated, if `false`, only the bonds in the crystallographic
  unit cell is calculated. Default is `true`.
  
`'conjugate'`
: Introduce the conjugate of the couplings (by exchanging the interacting
  `atom1` and `atom2`). Default is `false`.
  
### Output Arguments
  
`SS`
: structure with fields `iso`, `ani`, `dm`, `gen`, `bq`, `dip` and
  `all`. It describes the bonds between spins. Every field is a matrix,
              where every column is a coupling between two spins. The
              first 3 rows contain the unit cell translation vector
              between the interacting spins, the 4th and 5th rows contain
              the indices of the two interacting spins in the
              [spinw.matom](spinw_matom) list. The following rows contains the
              strength of the interaction. For isotropic exchange it is a
              single number, for DM interaction it is a column vector
              `[DMx; DMy; DMz]`, for anisotropic interaction `[Jxx; Jyy;
              Jzz]` and for general interaction `[Jxx; Jxy; Jxz; Jyx; Jyy;
              Jyz; Jzx; Jzy; Jzz]` and for biquadratic exchange it is also
              a single number. For example:
  ```matlab
  SS.iso = [dl_a; dl_b; dl_c; matom1; matom2; Jval]
  ```
  If `plotmode` is `true`, two additional rows are added to `SS.all`,
              that contains the `idx` indices of the
              `obj.matrix(:,:,idx)` corresponding matrix for each
              coupling and the `idx` values of the couplings (stored in
              `spinw.coupling.idx`). The `dip` field contains the dipolar
              interactions that are not added to the `SS.all` field.
 
`SI`
: single ion properties stored in a structure with fields:
  * `aniso`   Matrix with dimensions of $$[3\times 3\times n_{magAtom}]$$,
              where the classical energy of the $$i$$-th spin is expressed
              as `E_aniso = spin(:)*A(:,:,i)*spin(:)'`
	* `g`       g-tensor, with dimensions of $$[3\times 3\times n_{magAtom}]$$. It determines
              the energy of the magnetic moment in external field:
              `E_field = B(:)*g(:,:,i)*spin(:)'`
	* `field`   External magnetic field in a row vector with three elements $$(B_x, B_y, B_z)$$.
 
`RR`
: Positions of the atoms in lattice units in a matrix with dimensions of $$[3\times n_{magExt}]$$.
  
### See Also
  
[spinw.table](spinw_table) \| [spinw.symop](spinw_symop)
 
*[DM]: Dzyaloshinskii-Moriya
 

{% include links.html %}
