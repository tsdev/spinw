---
{title: sw_cff, link: sw_cff, summary: returns the atomic charge form factor values
    for X-ray scattering, keywords: sample, sidebar: sw_sidebar, permalink: sw_cff,
  folder: swfiles, mathjax: 'true'}

---
  
### Syntax
  
`[formfactval, coeff] = sw_cff(atomname, {q})`
  
### Description
  
`[formfactval, coeff] = sw_cff(atomname, {q})` returns the atomic charge
form factors for x-ray scattering. The provided form factor values at Q=0
are normalized to Z.
  
### Input Arguments
  
`atomName`
: String, contains the name of the ion in using the symbol of
  the element following the charge, e.g. `'Cr3+'`. It can be also
  the coefficients to calculate f. If the string contains
  whitespace, the first word will be used as input.
  
`Q`
: Momentum transfer in Å$$^{-1}$$ units in a matrix with dimensions of
  $$[1\times n_Q]$$ or $$[3\times n_Q]$$, optional.
  
### Output Arguments
  
`formFactVal`
: Value of form factor, evaluated at the $$Q$$ points if $$Q$$ is
              defined.
 
`coeff`
: Form factor coefficients according to the following
              formula:
                  $$f_0(Q_s) = c + \sum_{i=1}^5 a_i\cdot \exp(-b_i Q_s^2)$$,
              where $$Q_s = \frac{Q}{4*\pi}$$ and $$[a_1, b_1, a_2, b_2, ... c]$$ are the
              coefficients.
  
### See Also
  
[sw_mff](sw_mff)
 

{% include links.html %}
