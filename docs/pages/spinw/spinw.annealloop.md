---
{title: spinw.annealloop method, link: spinw.annealloop, summary: performs simulated
    annealing on the magnetic structure and measurements, keywords: sample, sidebar: sw_sidebar,
  permalink: spinw_annealloop.html, folder: spinw, mathjax: 'true'}

---

### Syntax

`stat = annealloop(obj,Name,Value)`

### Description

The function can deal only with single ion anisotropy and isotropic
exchange interactions in 1, 2 or 3 spin dimensions. General and DM
interactions are not supported yet!
 

### Input Arguments

`obj`
: [spinw](spinw.html) object.

### Name-Value Pair Arguments

`'spinDim'`
: Dimensionality of the magnetic moments.
      1   Ising spins
      2   XY spins
      3   Heisenberg spins [default]
  For Ising (spinDim=1) and XY (spinDim=2) models only isotropic
  exchange interaction and magnetic field can be used. For Ising
  the direction of the spins are along x-axis, for XY model the
  the xy-plane. Magnetic fields perpendicular to these directions
  are omitted.

`'func'`
: Function that changes the parameters in the spinw object in every
  loop. Default function is to change the temperature:
      @(obj,T)obj.temperature(T)
  The function takes two input: spinw objec and a parameter vector.

`'x'`
: Matrix of values of the loop parameter, with dimensions of
  [nPar nStep]. Default is 1. In the i-th loop the loop function
  is called as:
      func(obj,x(:,i));

`'random'`
: Random initial conditions before the first loop, if initial
  spin configuration is undefined (obj.mag_str.S is empty) the
  initial configuration is automaticly random independently of
  the value of random. Default is false.

`'nMC'`
: Number of Monte-Carlo steps per spin at each loop. Default is
  100.

`'nORel'`
: Number of over-relaxation steps after every Monte-Carlo
  steps. It rotates the spins around the direction of the local
  field by 180deg. It is reversible and microcanonical if the
  single ion anisotropy is zero. Default is 0.

`'nStat'`
: Number of cycles at the end of each loop to calculate
  statistical averages. Default is 100.

`'boundary'`
: Boundary conditions of the extended unit cell.
      'free'  Free, interactions between extedned unit cells are
              omitted.
      'per'   Periodic, interactions between extended unit cells
              are retained.
  Default is {'per' 'per' 'per'}.

`'verbosity'`
: Controls output to the screen.
      0   suppresses all output
      1   gives final report only [default]
      2   plots temperature changes and final report

`'nExt'`
: The size of the magnetic cell in number of unit cells, to
  provide input information to 'fStat'. Default is from
  obj.mag_str.N_ext.

`'fStat'`
: Function handle to evaluate after at the end of the
  cooling scedule during the last nStat Monte-Carlo steps. The
  function returns a single structure and takes fixed input
  parameters:
      struct = fStat(state, struct, T, E, M, nExt).
  The function is called once before the annealing process when
  state=1 to initialise the parameters. The function is called
  after every Monte-Carlo steps with state=2 and the output of
  the previous function call is assigned to the input struct.
  fStat is called once again in the end with state=3 to calculate
  final parameters (in the last run, input struct.param contains
  all the annealing parameters).
  Default is <a href="matlab: doc sw_fstat">@sw_fstat</a>.

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

`'saveObj'`
: If true, the spinw object is saved after every annealing step for
  debugging purposes. Default is false.

`'title'`
: Gives a title string to the simulation that is saved in the
  output.

### Output Arguments

stat      A struct type data that contains the calculated thermodynamical
          averages and the parameters of the simulation for evry value of
          X with the following fields:
param     All input parameter values of the annealloop function.
obj       The copy of the input spinw class obj with the final magnetic
          structure.
M         Components of the magnetisation after the last annealing
          run, dimensions are [3 nMagExt].
E         Magnetic energy of the system after the last annealing run.
T         Final temperature of the sample.
Depending on the 'fStat' parameter, additional fields are included. Using
the default function (@sw_fstat) the following parameters are calculated:
avgM      Average components of the magnetisation over nStat runs,
          dimensions are [3 nMagExt].
stdM      Standard deviation of the mgnetisation components over
          nStat runs, dimensions are [3 nMagExt].
avgE      Average system energy per spin over nStat runs, scalar.
stdE      Standard deviation of the system energy per spin over
          nStat runs, scalar.
Cp        Heat capacity of the sample: (<E^2>-<E>^2)/kB/T^2.
Chi       Magnetic susceptibility of the sample: (<M^2>-<M>^2)/kB/T.
 Reference:
   Kirkpatrick, S., Gelatt, C.D., & Vecchi, M.P. (1983). Optimization by
   Simulated Annealing. _Science, 220_, 671-680.

### See Also

[spinw](spinw.html) \| [spinw.optmagstr](spinw_optmagstr.html) \| [sw_fsub](sw_fsub.html) \| [sw_fstat](sw_fstat.html)

