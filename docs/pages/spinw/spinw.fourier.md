---
{title: spinw.fourier method, link: spinw.fourier, summary: calculates the Fourier
    transformation of the Hamiltonian, keywords: sample, sidebar: sw_sidebar, permalink: spinw_fourier,
  folder: spinw, mathjax: 'true'}

---
  
### Syntax
  
`F = fourier(obj,Q,Name,Value)`
  
### Description
  
`F = fourier(obj,hkl,Name,Value)` calculates the following Fourier sum:
 
$$J(\mathbf{k}) = \sum_{i,j} J_{i,j} * \exp(i \mathbf{k}\cdot \mathbf{d}_{i,j})$$
 
The code is optimised for calculating the sum for large number of wave
vectors and alternatively for a large number of $$d_{i,j}$$ vectors (large
system size). The single ion anisotropy is not included in the sum.
  
### Input Arguments
  
`obj`
: [spinw](spinw) object.
  
`Q`
: Defines the $$Q$$ points where the spectra is calculated, in reciprocal
  lattice units, size is $$[3\times n_{Q}]$$. $$Q$$ can be also defined by
  several linear scan in reciprocal space. In this case `Q` is cell type,
  where each element of the cell defines a point in $$Q$$ space. Linear scans
  are assumed between consecutive points. Also the number of $$Q$$ points can
  be specified as a last element, it is 100 by defaults. 
    
  For example to define a scan along $$(h,0,0)$$ from $$h=0$$ to $$h=1$$ using
  200 $$Q$$ points the following input should be used:
  ```matlab
  Q = {[0 0 0] [1 0 0]  50}
  ```
 
  For symbolic calculation at a general reciprocal space point use `sym`
  type input. 
 
  For example to calculate the spectrum along $$(h,0,0)$$ use:
  ```matlab
  Q = [sym('h') 0 0]
  ```
  To calculate spectrum at a specific $$Q$$ point symbolically, e.g. at
  $$(0,1,0)$$ use:
  ```matlab
  Q = sym([0 1 0])
  ```
  
### Name-Value Pair Arguments
  
`'extend'`
: If `true`, the Fourier transform will be calculated on the
  magnetic supercell, if `false` the crystallographic cell will
  be considered. Default is `true.`
  
`'isomode'`
: Defines how Heisenberg/non-Heisenberg Hamiltonians are
  treated. Can have the following values:
  * `'off'`   Always output the $$[3\times 3]$$ form of the
              Hamiltonian, (default).
  * `'auto'`  If the Hamiltonian is Heisenberg, only output
              one of the diagonal values from the $$[3\times 3]$$
              matrices to reduce memory consumption.
  
`'fid'`
: Defines whether to provide text output. Default is defined
  by the `swpref.getpref('fid')` command. The possible values:
  * `0`       No text output is generated.
  * `1`       Text output in the MATLAB Command Window.
  * `fid`     File ID provided by the [fopen](https://www.mathworks.com/help/matlab/ref/fopen.html) command, the
              output is written into the opened file stream.
  
### Output Arguments
  
`res` struct type with the following fields:
* `ft`        contains the Fourier transform in a matrix with dimensions
              $$[3\times 3\times n_{magExt}\times n_{magExt}\times
              n_{hkl}]$$ or $$[1\times 1\times n_{magExt}\times n_{magExt}\times n_{hkl}]$$
              for Heisenberg and non-Heisenberg Hamiltonians respectively
              (if isomode is `'auto'`). Here $$n_{magExt}$$ is the number of
              magnetic atoms in the magnetic cell and $$n_{hkl}$$ is the number
              of reciprocal space points.
* `hkl`       Matrix with the given reciprocal space points stored in a
              matrix with dimensions $$[3\times n_{hkl}]$$.
* `isiso`     True is the output is in Heisenberg mode, when the `ft`
              matrix has dimensions of $$[1\times 1\times n_{magExt}\times n_{magExt}\times n_{hkl}]$$,
              otherwise it is `false`.
  
### See Also
  
[spinw.optmagk](spinw_optmagk)
 

{% include links.html %}
