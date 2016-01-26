function fitpow(obj,varargin)
% fits spin wave spectra to experimental spectral data
%
% fitsp = FITPOW(obj, 'Option1', Value1, ...)
%
% Options:
%
% func      Function to change the Hamiltonian in obj, it has the following
%           header:
%                    obj = @func(obj, x);
% datapath  Path to the file that stores the experimental data. For the
%           input data format see <a href="matlab:doc sw_readspec">sw_readspec</a>.
% Evect     Vector, defines the energy binning of the calculated
%           dispersion. Larger binning steps solve the issue of fitting
%           unresolved modes. Size is [1 nE].
% xmin      Minimum limit of the optimisation parameters, optional.
% xmax      Maximum limit of the optimisation parameters, optional.
% x0        Starting value of the optimisation parameters. If empty
%           or undefined, then random values are used.
% nRun      Number of consecutive fitting runs, each result is saved in the
%           output fitsp.x and R arrays. If the Hamiltonian given by the
%           random x parameters is incompatible with the ground state,
%           those x values will be skipped and new random x values will be
%           generated. Default is 1.
% nMax      Maximum number of runs, including the ones that produce error
%           (due to incompatible ground state). Default is 1000.
% hermit    Method for matrix diagonalization:
%                  true      J.H.P. Colpa, Physica 93A (1978) 327,
%                  false     R.M. White, PR 139 (1965) A450.
%           Colpa: the grand dynamical matrix is converted into another
%                  Hermitian matrix, that will give the real eigenvalues.
%           White: the non-Hermitian g*H matrix will be diagonalised,
%                  that is not the elegant method.
%           Advise:
%           Always use Colpa's method that is faster, except when small
%           imaginary eigenvalues are expected. In this case only White's
%           method work.
%           Default is true.
% epsilon   Small number that controls wether the magnetic structure is
%           incommensurate or commensurate, default value is 1e-5.
%
% Parameters for visualizing the fit results:
%
% plot      If true, the measured dispersion is plotted together with the
%           fit. Default is true.
% iFact     Factor of the plotted simulated spin wave intensity (red
%           ellipsoids).
% lShift   Vertical shift of the Q point labels on the plot.
%
%
% Optimisation options:
%
% tolx          Minimum change of x when convergence reached, default
%               value is 1e-4.
% tolfun        Minimum change of the R value when convergence reached,
%               default value is 1e-5.
% maxfunevals   Maximum number of function evaluations, default value is
%               1e7.
%
% Output:
%
% Output fitsp is struct type with the following fields:
% obj       Copy of the input sw class object, with the best fitted
%           Hamiltonian.
% x         Final values of the fitted parameters, dimensions are
%           [nRun nPar]. The rows of x are sorted according to increasing R
%           values.
% R         R-value, goodness of the fit, dimensions are [nRun 1], sorted
%           in increasing order.
% exitflag  Exit flag of the fminsearch command.
% output    Output of the fminsearch command.
%
% Any option used by SW.SPINWAVE function are also accepted.
%
% See also SW.SPINWAVE, SW_EGRID, SW_NEUTRON, SW_READSPEC, FMINSEARCH.
%


end