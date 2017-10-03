---
{title: spinw.fitspec method, link: spinw.fitspec, summary: fits spin wave spectra
    to experimental spectral data, keywords: sample, sidebar: sw_sidebar, permalink: spinw_fitspec,
  folder: spinw, mathjax: 'true'}

---

### Syntax

`fitsp = fitspec(obj,Name,Value)`

### Description

The function uses a heuristic method to fit spin wave spectrum using a
few simple rules to define the R-value of the fit:
  1 All calculated spin wave modes that are outside of the measured
    energy range will be omitted.
  2 Spin wave modes that are closer to each other than the given energy
    bin will be binned together and considered as one mode in the fit.
  3 If the number of calculated spin wave modes after applying rule 1&2 
    is larger than the observed number, the weakes simulated modes will
    be removed from the fit.
  4 If the number of observed spin wave modes is larger than the observed
    number, fake spin wave modes are added with energy equal to the
    limits of the scan; at the upper or lower limit depending on which is
    closer to the observed spin wave mode.
After these rules the number of observed and simulated spin wave modes
will be equal. The R-value is defined as:
 
      R = sqrt( nE^(-1) * sum_i_q (E_i_q_sim - E_i_q_meas)^2/sigma_i_q^2 ),
 
where (i,q) indexing the spin wave mode and momentum. E_sim and E_meas
are the simulated and measured spin wave energies, sigma is the standard
deviation of the measured spin wave energy determined previously by
fitting the inelastic peak. nE is the number of energies to fit.
 
The R value is optimized using particle swarm algorithm to find the
global minimum.
 

### Name-Value Pair Arguments

`'func'`
:Function to change the Hamiltonian in obj, it has the following
 header:
          obj = @func(obj, x);

`'datapath'`
:Path to the file that stores the experimental data. For the
 input data format see <a href="matlab:doc sw_readspec">sw_readspec</a>.

`'Evect'`
:Vector, defines the energy binning of the calculated
 dispersion. Larger binning steps solve the issue of fitting
 unresolved modes. Size is [1 nE].

`'xmin'`
:Minimum limit of the optimisation parameters, optional.

`'xmax'`
:Maximum limit of the optimisation parameters, optional.

`'x0'`
:Starting value of the optimisation parameters. If empty
 or undefined, then random values are used.

`'optimizer'`
:String that determines the optimizer to use, possible values:
     'pso'       Particle-swarm optimizer, see ndbase.pso,
                 default.
     'simplex'   Matlab built-in simplex optimizer, see
                 fminsearch.

`'nRun'`
:Number of consecutive fitting runs, each result is saved in the
 output fitsp.x and R arrays. If the Hamiltonian given by the
 random x parameters is incompatible with the ground state,
 those x values will be skipped and new random x values will be
 generated. Default is 1.

`'nMax'`
:Maximum number of runs, including the ones that produce error
 (due to incompatible ground state). Default is 1000.

`'hermit'`
:Method for matrix diagonalization:
        true      J.H.P. Colpa, Physica 93A (1978) 327,
        false     R.M. White, PR 139 (1965) A450.
 Colpa: the grand dynamical matrix is converted into another
        Hermitian matrix, that will give the real eigenvalues.
 White: the non-Hermitian g*H matrix will be diagonalised,
        that is not the elegant method.
 Advise:
 Always use Colpa's method that is faster, except when small
 imaginary eigenvalues are expected. In this case only White's
 method work.
 Default is true.

`'epsilon'`
:Small number that controls wether the magnetic structure is
 incommensurate or commensurate, default value is 1e-5.

`'imagChk'`
:Checks that the imaginary part of the spin wave dispersion is
 smaller than the energy bin size. Default is true.

`'Parameters'`
: for visualizing the fit results:

`'plot'`
:If true, the measured dispersion is plotted together with the
 fit. Default is true.

`'iFact'`
:Factor of the plotted simulated spin wave intensity (red
 ellipsoids).

`'lShift'`
:ertical shift of the Q point labels on the plot.

`'Optimisation'`
:on options:

`'TolX'`
:    Minimum change of x when convergence reached, default
     value is 1e-4.

`'TolFun'`
:    Minimum change of the R value when convergence reached,
     default value is 1e-5.

`'MaxFunEvals'`
:s   Maximum number of function evaluations, default value is
     1e7.

`'MaxIter'`
:    Maximum number of iterations for the ndbse.pso optimizer.
     Default value is 20.

### Output Arguments

Output fitsp is struct type with the following fields:
obj       Copy of the input spinw class object, with the best fitted
          Hamiltonian.
x         Final values of the fitted parameters, dimensions are
          [nRun nPar]. The rows of x are sorted according to increasing R
          values.
R         R-value, goodness of the fit, dimensions are [nRun 1], sorted
          in increasing order.
exitflag  Exit flag of the fminsearch command.
output    Output of the fminsearch command.
Any option used by SPINW.SPINWAVE function are also accepted.

### See Also

[spinw.spinwave](spinw_spinwave) \| [spinw.matparser](spinw_matparser) \| [sw_egrid](sw_egrid) \| [sw_neutron](sw_neutron) \| [sw_readspec](sw_readspec)

{% include links.html %}
