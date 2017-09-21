---
{title: spinw.optmagsteep method, link: spinw.optmagsteep, summary: optimise magnetic
    structure using the method of steepest descent, keywords: sample, sidebar: sw_sidebar,
  permalink: spinw_optmagsteep.html, folder: spinw, mathjax: 'true'}

---

### Syntax

`optm = optmagsteep(obj,Name,Value)`

### Description

 

### Input Arguments

`obj`
: [spinw](spinw.html) object.

### Name-Value Pair Arguments

`'nRun'`
: Number of iterations, default value is 100 (it is usually enough).

`'boundary'`
: Boundary conditions of the extended unit cell.
      'free'  Free, interactions between extedned unit cells are
              omitted.
      'per'   Periodic, interactions between extended unit cells
              are retained.
  Default is {'per' 'per' 'per'}.

`'nExt'`
: The size of the magnetic cell in number of unit cells, to
  provide input information to 'fStat'.
  Default is from obj.mag_str.N_ext.

`'fSub'`
: Function to define sublattices for Monte-Carlo speedup.
  cGraph = fSub(conn,nExt), where cGraph is a (1,nMagExt) sized
  vector, conn is a (2,nConn) size matrix and nExt is equal to
  'nExt'. Default is <a href="matlab: doc sw_fsub">@sw_fsub</a>

`'subLat'`
: Vector that assigns all magnetic moments into non-interacting
  sublattices, contains a single index (1,2,3...) for every
  magnetic moment, size is (1,nMagExt). If undefined, the
  function defined in 'fSub' will be used to partition the
  lattice.

`'random'`
: Random initial conditions, if initial spin configuration
  is undefined (obj.mag_str.S is empty) the initial configuration
  is automaticly random independently of the value of random.
  Default is false.

`'TolX'`
: Minimum change of the magnetic moment when the algorithm stops.

`'saveAll'`
: Save moment directions for every loop, default value is false.

`'Hmin'`
: Minimum field value on the spin that moves the spin. If the
  molecular field absolute value is below this, the spin won't be
  turned. Default is zero.

`'plot'`
: If true, plot magnetic structure in real time. Default is false. 

`'pause'`
: Time in second to pause after every optimization loop to make
  slower movie. Default is 0.

### Output Arguments

'optm' is a struct type variable with the following fields:
obj       spinw object that contains the optimised magnetic structure.
M         Magnetic moment directions with dimensions [3 nMagExt], if
          'saveAll' parameter is true, it contains the magnetic structure
          after every loop in a matrix with dimensions [3 nMagExt nLoop].
dM     	The change of magnetic moment vector averaged over all moments
          in the last loop.
e         Energy per spin in the optimised structure.
param     Input parameters, stored in a struct.
nRun      Number of loops executed.
datestart Starting time of the function.
dateend   End time of the function.
title     Title of the simulation, given in the input.

### See Also

[spinw](spinw.html) \| [spinw.anneal](spinw_anneal.html) \| [sw_fsub](sw_fsub.html) \| [sw_fstat](sw_fstat.html)

{% include links.html %}
