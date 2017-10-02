---
{title: sw_tofres, link: sw_tofres, summary: includes Q resolution to the spectrum,
  keywords: sample, sidebar: sw_sidebar, permalink: sw_tofres, folder: swfiles, mathjax: 'true'}

---

### Syntax

`spectra = sw_tofres(spectra,Name,Value)`

### Description

Simulates the finite bin size of the cuts of TOF data.
 

### Input Arguments

`spectra`
: Input structure, contains calculated correlation functions
  withouth the resolution effect.

### Name-Value Pair Arguments

`'method'`
: String, determines the method to genera the Q points, options:
      'random'    The bin volume will be randomly sampled.
      'grid'      The bin volume will be split up to a regular
                  grid.

`'dQ'`
: Vector with three numbers of scalar. The width of the Q bin
  along the three reciprocal lattice directions. The spectrum
  will be integrated in the Q+/-(dQ/2) range. DEfault value is
  [0.1 0.1 0.1].

`'nQ'`
: Vector with three numbers or scalar. Gives the number of Q
  points along the three reciprocal lattice directions to average
  over or the number of random Q points for the random method.

### Output Arguments

spectra that contains the calculated intensity in the swConv field.

### See Also

[sw_egrid](sw_egrid) \| [sw_instrument](sw_instrument)

{% include links.html %}
