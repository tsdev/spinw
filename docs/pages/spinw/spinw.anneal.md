---
{title: spinw.anneal method, link: spinw.anneal, summary: performs simulated annealing
    of spins, keywords: sample, sidebar: sw_sidebar, permalink: spinw_anneal, folder: spinw,
  mathjax: 'true'}

---
  
### Syntax
  
`stat = anneal(obj,Name,Value)`
  
### Description
  
`stat = anneal(obj,Name,Value)` performs simulated annealing on the spin
Hamiltonian defined in `obj`. It assumes a classical spin system and
employs the [Metropolis–Hastings
algorithm](https://en.wikipedia.org/wiki/Metropolis–Hastings_algorithm)
for state updates. The annealing is performed from a given initial
temperature down to final temperature with user defined steps and number
of Monte-Carlo cycles per temperature. The `spinw.anneal` can also
measure different thermodynamic quantities, such as heat capacity. The
function can deal with single ion anisotropy and arbitrary exchange
interactions. It can also restrict the spin degrees of freedom from 3
(default) to 2 or 1 to speed up simulation on xy and Ising systems. For
these restricted dimensions only isotropic exchange interactions are
allowed. Also the g-tensor is assumed to be 2.
   
{% include warning.html content=" The calculated energies does not contain the self energy (spin
coupled to itself), thus the energies calculated here can differ from the
output of [spinw.energy](spinw_energy)." %}
   
  
### Input Arguments
  
`obj`
: [spinw](spinw) object.
  
### Name-Value Pair Arguments
  
`'spinDim'`
: Dimensionality of the spins:
  * **1**     Ising spins.
  * **2**     XY spins.
  * **3**     Heisenberg spins, default.
 
  For Ising (spinDim=1) and XY (spinDim=2) models only isotropic
  exchange interaction and magnetic field can be used. For Ising
  the direction of the spins are along x-axis, for XY model the
  the xy-plane. Magnetic fields perpendicular to these directions
  are omitted.
  
`'initT'`
: The initial temperature, can be any positive number
  in units of Kelvin. Default value is 1.
  
`'endT'`
: Temperature at which the annealing will stop, can be any positive number
  smaller than `initT`, unit is Kelvin.
  Default value is $$10^{-3}$$.
  
`'cool'`
: Defines how the following temperature value is calculated from the
  previous one using a user function. Any function that takes a scalar as input and
  returns a smaller but positive scalar as output. Default is `@(T)(0.92*T)`.
  
`'random'`
: If `true` the initial spin orientation will be random, which is
  effectively a $$T=\infty$$ paramagnet. If the initial spin configuration
  is undefined (`obj.magstr.S` is empty) the initial configuration
  is always random independently of this parameter.
  Default value is `false`.
  
`'nMC'`
: Number of Monte-Carlo steps per spin at each temperature
  step  which includes thermalization and the sampling for the extracted 
  TD properties at the last temperature). Default is 100.
  
`'nORel'`
: Number of over-relaxation steps after every Monte-Carlo
  steps. It rotates the spins around the direction of the local field by
  180°. It is reversible and microcanonical if the single ion
  anisotropy is 0. Default value is 0, to omit over-relaxation.
  
`'nStat'`
: Number of cycles at the last temperature to calculate
  statistical averages. It has to be smaller or equal $$n_{MC}$$.
  Default value is 100.
  
`'boundary'`
: Boundary conditions of the extended unit cell:
  * **free**  Free, interactions between extedned unit cells are
              omitted.
  * **per**   Periodic, interactions between extended unit cells
              are retained.
  Default value is `{'per' 'per' 'per'}`.
  
`'verbosity'`
: Controls the amount of output to the Command Window:
  * **0**   Suppresses all output.
  * **1**   Gives final report only.
  * **2**   Plots temperature changes and final report, default value.
  
`'nExt'`
: The size of the magnetic cell in number of unit cells that can override
  the value stored in `obj.magstr.N_ext`, given by a row vector with
  three integers
  
`'fStat'`
: Function handle to measure TD quantities at the final temperature
  for `nStat` Monte-Carlo steps. The function returns a single structure
  and takes fixed input parameters:
  ```matlab
  parOut = fStat(state, parIn, T, E, M, nExt).
  ```
  The function is called once before the annealing process
  when `state=1` to initialise the parameters. The function is called
  after every Monte-Carlo cycle with `state=2` and the output of the
  previous function call is assigned to the input struct. `fStat` is called
  once again in the end with `state=3` to calculate final parameters (in
  the last run, input `parIn.param` contains all the annealing
  parameters). For comparison see the defaul function [sw_fstat](sw_fstat).
  Default value is `@sw_fstat`.
  
`'fSub'`
: Function to define sublattices for Monte-Carlo speedup. Function handle
  with the following header:
  ```matlab
  cGraph = fSub(conn,nExt)
  ```
  where `cGraph` is a row vector with $$n_{magExt}$$ number of elements
  `conn` is a matrix with dimensions of $$[2\times n_{conn}]$$ $$n_{ext}$$ is
  equal to `nExt`. For the SpinW implementation see [sw_fsub](sw_fsub). Default
  value is `@sw_fsub`.
  
`'subLat'`
: Vector that assigns all magnetic moments into non-interacting
  sublattices, contains a single index $$(1,2,3...)$$ for every
  magnetic moment, row vector with $$n_{magExt}$$ number of elements. If
  undefined, the function defined in `fSub` parameter will be used to
  partition the lattice.
  
`'title'`
: Gives a title string to the simulation that is saved in the
  output.
  
`'autoK'`
: Bin length of the autocorrelation vector. Should be a few times
  smaller than `nMC`. Default value is 0 when no autocorrelation function
  is calculated.
  
### Output Arguments
  
`stat` struct that contains the calculated TD averages and the parameters
of the simulation with the following fields:
 
`param`
: All input parameter values of the anneal function.
  
`obj`
: The clone of the input `obj` updated with the final magnetic
  structure.
  
`M`
: Components of the magnetisation after the last annealing
  run stored in a matrix with dimensions of $$[3\times n_{magExt}]$$.
  
`E`
: Energy of the system after the last annealing run, excluding the self
  energy.
  
`T`
: Final temperature of the sample.
  
Depending on the `fStat` parameter, additional fields are included. Using
the default function [sw_fstat](sw_fstat) the following parameters are also
calculated:
  
`avgM`
: Average value of the magnetisation vector sampled over `nStat` runs,
  stored in a matrix with dimensions of $$[3\times n_{magExt}]$$.
  
`stdM`
: Standard deviation of the mgnetisation vector sampled over
  `nStat` runs, stored in a matrix with dimensions of $$[3\times
  n_{magExt}]$$.
  
`avgE`
: Average system energy per spin averaged over `nStat` runs, scalar.
  
`stdE`
: Standard deviation of the system energy per spin over
  `nStat` runs, scalar.
  
`Cp`
: Heat capacity of the sample, calculated using the formula $$(\langle E^2\rangle-\langle E\rangle^2)/k_B/T^2$$.
  
`Chi`
: Magnetic susceptibility of the sample calculated using the formula $$(\langle M^2\rangle-\langle M\rangle^2)/k_B/T$$.
 
 
### Reference
 
   Kirkpatrick, S., Gelatt, C.D., & Vecchi, M.P. (1983). Optimization by
   Simulated Annealing. _Science, 220_, 671-680.
  
### See Also
  
[spinw](spinw) \| [spinw.optmagstr](spinw_optmagstr) \| [sw_fsub](sw_fsub) \| [sw_fstat](sw_fstat)
 
*[TD]: Thermodynamic
 

{% include links.html %}
