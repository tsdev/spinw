---
{title: sw_tofres( ), summary: includes Q resolution to the spectrum, keywords: sample,
  sidebar: sw_sidebar, permalink: sw_tofres.html, folder: swfiles, mathjax: 'true'}

---
includes Q resolution to the spectrum
 
spectra = SW_TOFRES(spectra, 'Option1', Value1, ...)
 
Simulates the finite bin size of the cuts of TOF data.
 
Input:
 
spectra   Input structure, contains calculated correlation functions
          withouth the resolution effect.
 
Options:
 
method    String, determines the method to genera the Q points, options:
              'random'    The bin volume will be randomly sampled.
              'grid'      The bin volume will be split up to a regular
                          grid.
dQ        Vector with three numbers of scalar. The width of the Q bin
          along the three reciprocal lattice directions. The spectrum
          will be integrated in the Q+/-(dQ/2) range. DEfault value is
          [0.1 0.1 0.1].
nQ        Vector with three numbers or scalar. Gives the number of Q
          points along the three reciprocal lattice directions to average
          over or the number of random Q points for the random method.
 
 
Output:
 
spectra that contains the calculated intensity in the swConv field.
 
See also SW_EGRID, SW_INSTRUMENT.
 
