---
{title: sw_econtract( ), link: sw_econtract, summary: 'converts (Q,omega) values to
    Qm values for diffraction instrument', keywords: sample, sidebar: sw_sidebar,
  permalink: sw_econtract.html, folder: swfiles, mathjax: 'true'}

---

### Syntax

`qm = sw_flatcone(q,Name,Value)`

### Description



### Input Arguments

`Q`
: Input values in reciprocal space in the scattering plane in
  Angstrom^-1 units, dimensions are [2 nQ].

### Name-Value Pair Arguments

`omega`
: Energy transfer value in meV, default is zero.

`lambda`
: Wavelength of the incident beam in Angstrom.

`ki`
: Momentum of the incidend beam in Angstrom^-1.

`sense`
: Scattering sense:
      1       detectors are on the right hand side from the
              incident beam direction, default.
     -1       detectors are on the left hand side from the
              incident beam direction.

### See Also

[sw_converter](sw_converter.html)

