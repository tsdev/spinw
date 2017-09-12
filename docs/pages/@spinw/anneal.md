---
{title: '@spinw.anneal( )', summary: performs simulated annealing on the magnetic
    structure, keywords: sample, sidebar: sw_sidebar, permalink: spinw_anneal.html,
  folder: '@spinw', mathjax: 'true'}

---
performs simulated annealing on the magnetic structure
 
stat = ANNEAL(obj, 'option1', value1 ...)
 
The function can deal only with single ion anisotropy and isotropic
exchange interactions in 1, 2 or 3 spin dimensions.
General and antisymmetric exchage interactions are not supported yet!
Also the g-tensor is fixed to 2.
 
WARNING!
The calculated energies doesn't contain the self energy (moment coupled
to itself), thus the energies calculated here can differ from the
result of the spinw.energy() function.
 
Input:
 
obj             Input object contains structural data, spinw class.
 
Options:
 
spinDim   Dimensionality of the magnetic moments.
              1   Ising spins
              2   XY spins
              3   Heisenberg spins [default]
          For Ising (spinDim=1) and XY (spinDim=2) models only isotropic
          exchange interaction and magnetic field can be used. For Ising
          the direction of the spins are along x-axis, for XY model the
          the xy-plane. Magnetic fields perpendicular to these directions
          are omitted.
initT     The initial temperature, can be any positive number,
          unit is Kelvin. Default is 1.
endT      Temperature at which to stop, can be any positive number
          smaller than 'InitTemp', unit is Kelvin.
          Default is 1e-3.
cool      Generates a new temperature from the previous one.
          Any function handle that takes a scalar as input and
          returns a smaller but positive scalar as output.
          Default is @(T) (.92*T).
random    Random initial conditions, if initial spin configuration
          is undefined (obj.mag_str.S is empty) the initial configuration
          is automaticly random independently of the value of random.
          Default is false.
nMC       Number of Monte-Carlo steps per spin at each temperature
          step to reach thermal equilibrium. Default is 100.
nORel     Number of over-relaxation steps after every Monte-Carlo
          steps. It rotates the spins around the direction of the local
          field by 180deg. It is reversible and microcanonical if the
          single ion anisotropy is zero. Default is 0.
nStat     Number of cycles at the last temperature to calculate
          statistical averages. It has to be smaller or equal nMC.
          Default is 100.
boundary  Boundary conditions of the extended unit cell.
              'free'  Free, interactions between extedned unit cells are
                      omitted.
              'per'   Periodic, interactions between extended unit cells
                      are retained.
          Default is {'per' 'per' 'per'}.
verbosity Controls output to the screen.
              0   suppresses all output
              1   gives final report only
              2   plots temperature changes and final report [default]
nExt      The size of the magnetic cell in number of unit cells, to
          provide input information to 'fStat'.
          Default is from obj.mag_str.N_ext.
fStat     Function handle to evaluate after at the end of the
          cooling scedule during the last nStat Monte-Carlo steps.
          The function returns a single structure and takes fixed
          input parameters:
              struct = fStat(state, struct, T, E, M, nExt).
          The function is called once before the annealing process
          when state=1 to initialise the parameters. The function
          is called after every Monte-Carlo steps with state=2 and
          the output of the previous function call is assigned to
          the input struct. fStat is called once again in the end
          with state=3 to calculate final parameters (in the last
          run, input struct.param contains all the annealing
          parameters).
          Default is <a href="matlab: doc sw_fstat">@sw_fstat</a>.
fSub      Function to define sublattices for Monte-Carlo speedup.
          cGraph = fSub(conn,nExt), where cGraph is a (1,nMagExt) sized
          vector, conn is a (2,nConn) size matrix and nExt is equal to
          'nExt'. Default is <a href="matlab: doc sw_fsub">@sw_fsub</a>
subLat    Vector that assigns all magnetic moments into non-interacting
          sublattices, contains a single index (1,2,3...) for every
          magnetic moment, size is (1,nMagExt). If undefined, the
          function defined in 'fSub' will be used to partition the
          lattice.
title     Gives a title string to the simulation that is saved in the
          output.
autoK     Bin length of the autocorrelation vector. Should be a few times
          smaller than nMC. Default is zero, no autocorrelation function
          is calculated.
 
Output:
 
stat      Struct that contains the calculated thermodynamical
          averages and the parameters of the simulation with the
          following fields:
 
param     All input parameter values of the anneal function.
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
 
See also SPINW, SPINW.OPTMAGSTR, SW_FSUB, SW_FSTAT.
 
