---
{title: spinw.optmagsteep method, link: spinw.optmagsteep, summary: quench optimization
    of magnetic structure, keywords: sample, sidebar: sw_sidebar, permalink: spinw_optmagsteep,
  folder: spinw, mathjax: 'true'}

---
  
### Syntax
  
`optm = optmagsteep(obj,Name,Value)`
  
### Description
  
`optm = optmagsteep(obj,Name,Value)` determines the lowest energy
magnetic configuration within a given magnetic supercell and previously
fixed propagation (and normal) vector (see [spinw.optmagk](spinw_optmagk)). It
iteratively rotates each spin towards the local magnetic field thus
achieving local energy minimum. Albeit not guaranteed this method often
finds the global energy minimum. The methods works best for small
magnetic cells and non-frustrated structures. Its execution is roughly
equivalent to a thermal quenching from the paramagnetic state.
  
### Input Arguments
  
`obj`
: [spinw](spinw) object.
  
### Name-Value Pair Arguments
  
`'nRun'`
: Number of iterations, default value is 100 (it is usually enough). Each
  spin will be quenched `nRun` times or until convergence is reached.
  
`'boundary'`
: Boundary conditions of the magnetic cell, string with allowed values:
  * `'free'`  Free, interactions between extedned unit cells are
              omitted.
  * `'per'`   Periodic, interactions between extended unit cells
              are retained.
 
  Default value is `{'per' 'per' 'per'}`.
  
`'nExt'`
: The size of the magnetic cell in number of crystal unit cells.
  Default value is taken from `obj.mag_str.nExt`.
  
`'fSub'`
: Function that defines non-interacting sublattices for parallelization.
  It has the following header:
      `cGraph = fSub(conn,nExt)`, where `cGraph` is a row vector with
      $$n_{magExt}$$ number of elements,
  `conn` is a matrix with dimensions of $$[2\times n_{conn}]$$ size matrix and $$n_{ext}$$ is equal to
  the `nExt` parameter. Default value is `@sw_fsub`.
  
`'subLat'`
: Vector that assigns all magnetic moments into non-interacting
  sublattices, contains a single index $$(1,2,3...)$$ for every magnetic
  moment in a row vector with $$n_{magExt}$$ number of elements. If
  undefined, the function defined in `fSub` will be used to partition the
  lattice.
  
`'random'`
: If `true` random initial spin orientations will be used (paramagnet),
  if initial spin configuration is undefined (`obj.mag_str.F` is empty)
  the initial configuration will be always random. Default value is
  `false`.
  
`'TolX'`
: Minimum change of the magnetic moment necessary to reach convergence.
  
`'saveAll'`
: Save moment directions for every loop, default value is `false`.
  
`'Hmin'`
: Minimum field value on the spin that moves the spin. If the
  molecular field absolute value is below this, the spin won't be
  turned. Default is 0.
  
`'plot'`
: If true, the magnetic structure in plotted in real time. Default value
  is `false`.
  
`'pause'`
: Time in second to pause after every optimization loop to slow down plot
  movie. Default value is 0.
  
### Output Arguments
  
`optm`
: Struct type variable with the following fields:
  * `obj`         spinw object that contains the optimised magnetic structure.
  * `M`           Magnetic moment directions with dimensions $$[3\times n_{magExt}]$$, if
                  `saveAll` parameter is `true`, it contains the magnetic structure
                  after every loop in a matrix with dimensions $$[3\times n{magExt}\times n_{loop}]$$.
  * `dM`          The change of magnetic moment vector averaged over all moments
                  in the last loop.
  * `e`           Energy per spin in the optimised structure.
  * `param`       Input parameters, stored in a struct.
  * `nRun`        Number of loops executed.
  * `datestart`   Starting time of the function.
  * `dateend`     End time of the function.
  * `title`       Title of the simulation, given in the input.
  
### See Also
  
[spinw](spinw) \| [spinw.anneal](spinw_anneal) \| [sw_fsub](sw_fsub) \| [sw_fstat](sw_fstat)
 

{% include links.html %}
