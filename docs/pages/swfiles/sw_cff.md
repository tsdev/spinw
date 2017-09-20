---
{title: sw_cff( ), link: sw_cff, summary: returns the atomic charge form factor values
    for X-ray scattering, keywords: sample, sidebar: sw_sidebar, permalink: sw_cff.html,
  folder: swfiles, mathjax: 'true'}

---

### Syntax

`[formfactval, coeff] = sw_cff(atomname, {q})`

### Description

The provided form factor values at Q=0 are normalized to Z.
 

### Input Arguments

`atomName`
: String, contains the name of the ion in using the symbol of
  the element following the charge, e.g. Cr3+. It can be also
  the coefficients to calculate f. If the string contains
  whitespace, the first word will be used as input.

`Q`
: Momentum transfer in Å$$^{-1}$$ units with dimensions of
  [1 nQ] or [3 nQ], optional.

### Output Arguments

formFactVal   Value of form factor, evaluated at the Q points if Q is
              defined.
coeff         Form factor coefficients according to the following
              formula:
                  f0(Qs) = c + SUM a_i*EXP(-b_i*(Qs^2))
                               i=1,5
              where Qs = Q/(4*pi) and [a_1 b_1 a_2 b_2 ... c] are the
              coefficients.

### See Also

[sw_mff](sw_mff.html)

