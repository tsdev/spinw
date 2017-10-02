---
{title: sw_mff, link: sw_mff, summary: returns the magnetic form factor values and
    the coefficients, keywords: sample, sidebar: sw_sidebar, permalink: sw_mff, folder: swfiles,
  mathjax: 'true'}

---

### Syntax

`[formfactval, coeff, s] = sw_mff(atomname, {q})`

### Description



### Input Arguments

`atomName`
: String, contains the name of the magnetic ion in FullProf
  notation (e.g. Cr^3+ --> 'MCR3' or 'Cr3'). It can be also a
  vector of the 7 coefficients, see below. If the string
  contains whitespace, the first word will be used as input.
  Can be also a cell of strings to calculate multiple ions.

`Q`
: Momentum transfer in Å$$^{-1}$$ units with dimensions of
  [1 nQ] or [3 nQ], optional.

### Output Arguments

formFactVal   Value of form factor, evaluated at the Q points if Q is
              defined.
coeff         Form factor coefficients according to the following
              formula:
              <j0(Qs)> = A*exp(-a*Qs^2) + B*exp(-b*Qs^2) + C*exp(-c*Qs^2) + D*exp(-d*Qs^2) + E,
              where Qs = Q/(4*pi) and A, a, B, ... are the coefficients.
              The (D,d) coefficients can be zero.
S             Value of the spin quantum number ('spin' column in magion.dat).
The source for the form factor data are:
[1] A.-J. Dianoux and G. Lander, Neutron Data Booklet (2003).
[2] K. Kobayashi, T. Nagao, and M. Ito, Acta Crystallogr. A. 67, 473 (2011).

{% include links.html %}
