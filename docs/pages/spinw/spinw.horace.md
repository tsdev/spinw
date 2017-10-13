---
{title: spinw.horace method, link: spinw.horace, summary: spin wave calculator with
    interface to Horace, keywords: sample, sidebar: sw_sidebar, permalink: spinw_horace,
  folder: spinw, mathjax: 'true'}

---
  
### Syntax
  
`[w, s] = horace(obj, qh, qk, ql,Name,Value)`
  
### Description
  
`[w, s] = horace(obj, qh, qk, ql,Name,Value)` produces spin wave
dispersion and intensity for [Horace](http://horace.isis.rl.ac.uk).
  
### Examples
  
This example creates a `d3d` object, a square in $$(h,k,0)$$ plane and in
energy between 0 and 10 meV. Then calculates the inelastice neutron
scattering intensity of the square lattice antiferromagnet stored in
`cryst` and plots a cut between 4 and 5 meV using the Horace `plot`
command.
```matlab
cryst = sw_model('squareAF',1)
d3dobj = d3d(cryst.abc,[1 0 0 0],[0,0.02,2],[0 1 0 0],[0,0.02,2],[0 0 0 1],[0,0.1,10])
d3dobj = disp2sqw_eval(d3dobj,@cryst.horace,{'component','Sperp'},0.1)
plot(cut(d3dobj,[],[],[4 5]))
```
 
{% include image.html file="generated/spinw_h_1.png" alt="colorslider('delete')" %}
  
### Input Arguments
  
`obj`
: [spinw](spinw) object.
  
`qh`, `qk`, `ql`
: Reciprocal lattice vectors in reciprocal lattice units.
  
### Name-Value Pair Arguments
  
`'component'`
: Selects the previously calculated intensity component to be
  convoluted. The possible options are:
  * `'Sperp'` convolutes the magnetic neutron scattering
              intensity ($$\langle S_\perp \cdot S_\perp\rangle$$ expectation value).
              Default value.
  * `'Sab'`   convolutes the selected components of the spin-spin
              correlation function.
  For details see [sw_egrid](sw_egrid).
  
`'norm'`
: If `true` the spin wave intensity is normalized to mbarn/meV/(unit
  cell) units. Default is `false`.
  
`'dE'`
: Energy bin size, for intensity normalization. Use 1 for no
  division by `dE` in the intensity.
  
`'param'`
: Input parameters (can be used also within Tobyfit). Additional
  parameters (`'mat'`,`'selector'`) might be necessary, for details see
  [spinw.matparser](spinw_matparser). All extra parameters of `spinw.horace`
  will be forwarded to the [spinw.matparser](spinw_matparser) function before
  calculating the spin wave spectrum (or any user written parser
  function). For user written functions defined with the
  following header:
  ```matlab
  func(obj,param)
  ```
  the value of the param option will be forwarded. For user
  functions with variable number of arguments, all input options
  of `spinw.horace` will be forwarded. In this case it is recommended
  to use [sw_readparam](sw_readparam) function to handle the variable number
  arguments within `func()`.
  
`'parfunc'`
: Parser function of the `param` input. Default value is
  `@spinw.matparser` which can be used directly by Tobyfit. For user
  defined functions the minimum header has to be:
  ```matlab
  func(obj,param)
  ```
  where obj is an spinw type object, param is the parameter
  values forwarded from` spinw.horace` directly.
  
`'func'`
: User function that will be called after the parameters set on
  the [spinw](spinw) object. It can be used to optimize magnetic
  structure for the new parameters, etc. The input should be a
  function handle of a function with a header:
  ```matlab
  fun(obj)
  ```
  
### Output Arguments
  
`w`
: Cell that contains the spin wave energies. Every cell elements
          contains a vector of spin wave energies for the corresponding
          input $$Q$$ vector.
 
`s`
: Cell that contains the calculated element of the spin-spin
          correlation function. Every cell element contains a vector of
          intensities in the same order as the spin wave energies in `w`.
  
### See Also
  
[spinw](spinw) \| [spinw.spinwave](spinw_spinwave) \| [spinw.matparser](spinw_matparser) \| [sw_readparam](sw_readparam)
 

{% include links.html %}
