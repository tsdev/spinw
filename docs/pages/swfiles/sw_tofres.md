---
{title: sw_tofres, link: sw_tofres, summary: convolutes the spectrum with a Q bin,
  keywords: sample, sidebar: sw_sidebar, permalink: sw_tofres, folder: swfiles, mathjax: 'true'}

---
  
### Syntax
  
`spectra = sw_tofres(spectra,Name,Value)`
  
### Description
  
`spectra = sw_tofres(spectra,Name,Value)` simulates the finite bin size
of the cuts of direct TOF neutron scattering data. It calculates the
spectrum at multiple points within the given bin volume and sums them up.
The function is usefull if relatively large bins were used to analyse the
data due to low signal to noise ratio of the measurement.
  
### Input Arguments
  
`spectra`
: Input structure, contains calculated correlation functions
  withouth the resolution effect.
  
### Name-Value Pair Arguments
  
`'method'`
: String that determines the method to generate the $$Q$$ points, options:
  * `'random'`    The bin volume will be randomly sampled.
  * `'grid'`      The bin volume will be split up to a finer regular
                  grid.
  
`'dQ'`
: Row vector with 3 elements. The width of the $$Q$$ bin
  along the three reciprocal lattice directions. The spectrum
  will be integrated in the $$Q\pm (\delta Q/2)$$ range. Default value is
  `[0.1 0.1 0.1]`.
  
`'nQ'`
: Row vector with 3 elements when `method` is `grid` and gives the
  number of $$Q$$ points along the three reciprocal lattice directions to
  average over. When `method` is `random` it is a scalar that determines
  the number of random $$Q$$ points.
 
`'fid'`
: Defines whether to provide text output. Default value is defined in
  `obj.fid`. The possible values are: 
  * `0`   No text output is generated.
  * `1`   Text output in the MATLAB Command Window.
  * `fid` File ID provided by the `fopen` command, the output is written
          into the opened file stream.
 
`'tid'`
: Determines if the elapsed and required time for the calculation is
  displayed. The default value is determined by the `tid` preference
  stored in [swpref]. The following values are allowed (for more details
  see [sw_timeit](sw_timeit)):
  * `0` No timing is executed.
  * `1` Display the timing in the Command Window.
  * `2` Show the timing in a separat pup-up window.
  
### Output Arguments
  
`spectra`
: Same as the input except that it contains the calculated intensity in
  the `swConv` field.
  
### See Also
  
[sw_egrid](sw_egrid) \| [sw_instrument](sw_instrument)
 
*[TOF]: Time Of Flight
 

{% include links.html %}
