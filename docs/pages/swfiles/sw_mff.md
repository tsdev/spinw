---
{title: sw_mff, link: sw_mff, summary: returns the magnetic form factor values and
    coefficients, keywords: sample, sidebar: sw_sidebar, permalink: sw_mff, folder: swfiles,
  mathjax: 'true'}

---
  
### Syntax
  
`[~, coeff, s] = sw_mff(atomname)`
  
`[formfactval, coeff, s] = sw_mff(atomname,Q)`
 
### Description
  
`[~, coeff, s] = sw_mff(atomname)` returns the magnetic form
factor coefficients for the magnetic atom identified by a string, e.g.
`'MCR3'`. The function reads the `magion.dat` file for the stored form
factor coefficients.
 
`[formfactval, coeff, s] = sw_mff(atomname,Q)` also calculates the form
factor values at the given $$Q$$ points (in Å$$^{-1}$$ units.
 
The source of the form factor data are:
1. A.-J. Dianoux and G. Lander, Neutron Data Booklet (2003).
2. K. Kobayashi, T. Nagao, and M. Ito, Acta Crystallogr. A. 67, 473 (2011).
 
### Input Arguments
  
`atomName`
: String, contains the name of the magnetic ion in FullProf
  notation (e.g. for Cr$$^{3+} use `'MCR3'` or `'Cr3'`). It can be also a
  vector of the 7 form factor coefficients. If the string contains
  whitespace, the first word will be used as input. Can be also a cell of
  strings to calculate coefficients for multiple ions.
  
`Q`
: Momentum transfer in Å$$^{-1}$$ units in a matrix with dimensions of
  $$[1\times n_Q]$$ or $$[3\times n_Q]$$.
  
### Output Arguments
  
`formFactVal`
: Value of the form factor, evaluated at the given $$Q$$ points.
 
`coeff`
: Form factor coefficients according to the following formula:
    
  $$\langle j_0(Q_s)\rangle = A\exp(-a\cdot Q_s^2) + B\exp(-b\cdot Q_s^2) + C\exp(-c\cdot Q_s^2) + D\exp(-d\cdot Q_s^2) + E$$
 
  where $$Q_s = \frac{Q}{4\pi}$$ and $$A$$, $$a$$, $$B$$, ... are the coefficients.
  The $$D$$ and $$d$$ coefficients can be zero.
  
`S`
: Value of the spin quantum number (read from the spin column in `magion.dat`).
 

{% include links.html %}
