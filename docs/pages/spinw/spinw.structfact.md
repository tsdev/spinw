---
{title: spinw.structfact method, link: spinw.structfact, summary: calculates magnetic
    and nuclear structure factor, keywords: sample, sidebar: sw_sidebar, permalink: spinw_structfact,
  folder: spinw, mathjax: 'true'}

---
  
### Syntax
  
`sFact   = structfact(obj, kGrid,Name,Value)`
 
`sfTable = structfact(obj, kGrid,Name,Value)`
 
### Description
  
`sFact   = structfact(obj, kGrid,Name,Value)` returns the calculated
structure factors in units of barn. Magnetic structures (FM, AFM and
helical) are checked against
[FullProf](https://www.ill.eu/sites/fullprof/). The structure factor
includes the site occupancy and Debye-Waller factors calculated from
`obj.unit_cell.biso`, using the same definition as in FullProf.
  
### Input Arguments
  
`obj`
: [spinw](spinw) object.
  
`kGrid`
: Defines the reciprocal lattice vectors where the structure
     factor is to be calculated. For commensurate structures these
     are the possible positions of the magnetic Bragg peaks. For
     incommensurate helical/conical structures 3 Bragg peaks
     positions are possible: $$(\mathbf{k}-\mathbf{k}_m,\mathbf{k},\mathbf{k}+\mathbf{k}_m) around every reciprocal
     lattice vector. In this case still the integer positions have
     to be given and the code calculates the intensities at all
     three points.
  
### Name-Value Pair Arguments
  
`'mode'`
: String, defines the type of calculation:
  * `mag`     Magnetic structure factor and intensities for
              unpolarised neutron scattering.
  * `nucn`    Nuclear structure factor and neutron scattering
              intensities.
  * `nucx`    X-ray scattering structure factor and
              intensities.
  
`'sortq'`
: Sorting the reflections according to increasing momentum
  value if `true`. Default is `false`.
  
`'formfact'`
: If true, the magnetic form factor is included in the structure factor
  calculation. The form factor coefficients are stored in
  `obj.unit_cell.ff(1,:,atomIndex)`. Default value is `false`.
 
`'formfactfun'`
: Function that calculates the magnetic form factor for given $$Q$$ value.
  value. Default value is `@sw_mff`, that uses a tabulated coefficients
  for the form factor calculation. For anisotropic form factors a user
  defined function can be written that has the following header:
  ```matlab
  F = formfactfun(atomLabel,Q)
  ```
  where the parameters are:
  * `F`           row vector containing the form factor for every input 
                  $$Q$$ value
  * `atomLabel`   string, label of the selected magnetic atom
  * `Q`           matrix with dimensions of $$[3\times n_Q]$$, where each
                  column contains a $$Q$$ vector in $$Å^{-1}$$ units.
 
`'gtensor'`
: If true, the g-tensor will be included in the structure factor
  calculation. Including anisotropic g-tensor or different
  g-tensor for different ions is only possible here.
 
`'lambda'`
: Wavelength. If given, the $$2\theta$$ value for each reflection
  is calculated.
  
`'dmin'`
: Minimum $$d$$-value of a reflection, all higher order
  reflections will be removed from the results.
  
`'output'`
: String, defines the type of the output:
  * `struct`  Results are returned in a struct type variable,
              default.
  * `table`   Results are returned in a table type output for
              easy viewing and exporting.
  
`'tol'`
: Tolerance of the incommensurability of the magnetic
  ordering wavevector. Deviations from integer values of the
  ordering wavevector smaller than the tolerance are considered
  to be commensurate. Default value is $$10^{-4}$$.
  
`'fitmode'`
: Speed up the calculation for fitting mode (omitting
  cloning the [spinw](spinw) object into the output). Default is `false`.
  
### Output Arguments
  
`sFact`
: Structure with the following fields:
   * `F2`     Magnetic structure factor in a matrix with dimensions
              $$[3\times n_{hkl}]$$.
   * `Mk`     Square of the 3 dimensional magnetic structure factor,
              dimensions are:
              $$[n_{ext}(1)\cdot f_{ext}(1)\times n_{ext}(2)\cdot f_{ext}(2)\times n_{ext}(3)\cdot f_{ext}(3)]$$,
              where $$n_{ext}$$ is the size of the extended unit cell.
   * `hkl`    Contains the input $$Q$$ values in a matrix with dimensins of $$[3\times n_{hkl}]$$.
   * `hklA`   Same as `hkl`, but in Å$$^{-1}$$ units in the
              $$xyz$$ Cartesian coordinate system.
   * `incomm` Whether the spectra calculated is incommensurate or not.
   * `formfact` Cell containing the labels of the magnetic ions if form
              factor in included in the spin-spin correlation function.
   * `{tth}`  $$2\theta$$ value of the reflection for the given wavelength,
              only given if a wavelength is provided.
   * `obj`    Clone of the input `obj` object.
 
`sfTable`
: Table, optional output for quick viewing and saving the output into a
  text file.
  
### See Also
  
[sw_qgrid](sw_qgrid) \| [sw_plotsf] \| [sw_intsf] \| [spinw.anneal](spinw_anneal) \| [spinw.genmagstr](spinw_genmagstr)
 

{% include links.html %}
