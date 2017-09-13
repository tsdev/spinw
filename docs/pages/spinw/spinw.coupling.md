---
{title: spinw.coupling( ), summary: Stores the list of spin-spin couplings., keywords: sample,
  sidebar: sw_sidebar, permalink: spinw_coupling.html, folder: spinw, mathjax: 'true'}

---
Stores the list of spin-spin couplings.
Sub fields are:
  dl      distance between the unit cells of two interacting
          spins, stored in a 3 x nCoupling matrix
  atom1   first magnetic atom, pointing to the list of
          magnetic atoms in spinw.matom list, stored in a
          1 x nCoupling vector
  atom2   second magnetic atom, stored in a  1 x nCoupling
          vector
  mat_idx stores pointers to matrices for every coupling in a
          3 x nCoupling matrix, maximum three matrices per
          coupling (zeros for no coupling)
  idx     increasing indices for the symmetry equivalent
          couplings, starting with 1,2,3...
  type    Type of coupling corresponding to mat_idx matrices.
          Default is 0 for quadratic exchange. type = 1 for
          biquadratic exchange.
  sym     If true the bond symmetry operators will be applied
          on the assigned matrix.
  rdip    Maximum distance until the dipolar interaction is
          calculated. Zero vlue means no dipolar interactions
          are considered.
  nsym    The largest bond 'idx' value that is generated
          using the space group operators. Typically very
          long bonds for dipolar interactions won't be
          calculated using space group symmetry.
 
See also SPINW.GENCOUPLING, SPINW.ADDCOUPLING, SPINW.FIELD.
