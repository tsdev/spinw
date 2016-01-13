function obj = horace_setpar(obj, varargin)
% Setup sw object to be called from Horace.
%
% HORACE_SETPAR(obj)
%
% Used with sw.horace to get spin wave dispersion and intensity for Horace
% (<a href=http://horace.isis.rl.ac.uk>http://horace.isis.rl.ac.uk</a>).
%
% Input:
%
% obj       Input sw object.
%
% Output:
%
% The function sets parameters in the 'horace' field of the obj sw object.
%
% Options:
%
% component Selects the previously calculated intensity component to be
%           convoluted. The possible options are:
%               'Sperp' convolutes the magnetic neutron scattering
%                       intensity (<Sperp * Sperp> expectation value).
%                       Default.
%               'Sab'   convolutes the selected components of the spin-spin
%                       correlation function. Letter a and b can be 'x',
%                       'y' or 'z'. For example: 'Sxx' will convolute the
%                       xx component of the correlation function with the
%                       dispersion. xyz is the standard coordinate system,
%                       see online documentation of SpinW.
%           Any linear combination of the above are allowed, for example:
%           'Sxx+2*Syy' convolutes the linear combination of the xx
%           component of the spin-spin correlation function and the yy
%           component.
% norm      If true the spin wave intensity is normalized to mbarn/meV/(unit
%           cell) units. Default is false.
% dE        Energy bin size, for intensity normalization. Use 1 for no
%           division by dE in the intensity.
% tol       Tolerance of the incommensurability of the magnetic
%           ordering wavevector. Deviations from integer values of the
%           ordering wavevector smaller than the tolerance are
%           considered to be commensurate. Default value is 1e-4.
% omega_tol Tolerance on the energy difference of degenerate modes when
%           diagonalising the quadratic form, default is 1e-5.
% optmem    Parameter to optimise memory usage. The list of hkl values
%           will be cut into optmem number of pieces and will be
%           calculated piece by piece to decrease memory usage. Default
%           of optmem is zero, when the number of slices are determined
%           automatically from the available free memory.
% hermit    Method for matrix diagonalization:
%                  true      J.H.P. Colpa, Physica 93A (1978) 327,
%                  false     R.M. White, PR 139 (1965) A450.
%           Colpa: the grand dynamical matrix is converted into another
%                  Hermitian matrix, that will give the real
%                  eigenvalues.
%           White: the non-Hermitian g*H matrix will be diagonalised,
%                  that is not the elegant method.
%           Advise:
%           Always use Colpa's method, except when small imaginary
%           eigenvalues are expected. In this case only White's method
%           work. The solution in this case is wrong, however by
%           examining the eigenvalues it can give a hint where the
%           problem is.
%           Default is true.
% notwin    If true, the spectra of the twins won't be calculated.
%           Default is false.
% param     Input vector P with nPar elements that contains the
%           new values to be assignd to elements of sw.matrix.mat
%           matrix.
% mapping   A mapping between values in the 'param' vector and matrices
%           in the sw object. To select matrices with given labels use a
%           cell of strings with nPar elements, for example
%           M = {'J1','J2'}. This will change the diagonal elements of
%           matrices J1 and J2 to a given value that is provided in the
%           'param' option. Alternatively the index of the matrices can
%           be given in a vector, such as [1 2] (index runs according
%           to the order of the previous creation of the matrices using
%           sw.addmatrix).
%           To assign parameter value only to a selected element of a
%           3x3 matrix, a bracket notation can be used in any string,
%           such as 'D(3,3)', in this case only the (3,3) element of
%           the 3x3 matrix of 'D' will be modified, the other elements
%           will be unchanged. To modify multiple elements of a matrix
%           at once, use the option 'selector'.
% selector  Matrix with dimensions of [3 3 nPar]. Each S(:,:,ii)
%           submatrix can contain +/-1 and 0. Where S(:,:,ii) contains
%           ones, the corresponding matrix elements of
%           sw.matrix.mat(:,:,M(ii)) will be changed to the value
%           P(ii)*S(:,:,ii) where P(ii) is the corresponding parameter
%           value. For example do assign a Dzyaloshinskii-Moriya vector
%           to the 'DM' matrix, the following input would be
%           sufficient:
%           P = [0.2 0.35 3.14]
%           M = {'DM' 'DM' 'DM'}
%           S = cat(3,[0 0 0;0 0 1;0 -1 0],[0 0 -1;0 0 0;1 0 0],[0 1 0;-1 0 0;0 0 0])
%           sw.horace_setpar('mapping',M,'selector',S)
% init      Initialize the matrices of sw.matrix.mat with zeros for all
%           selected labels before assigning paramter values. Default
%           is false.
% func      Parser function of the 'param' input. Default is
%           @sw.matparser which can be used directly by Tobyfit. For user
%           defined functions the minimum header has to be:
%               func(obj,param)
%           where obj is an sw type object, param is the parameter values
%           forwarded from sw.horace directly.
% formfact  Setting, that determines whether the magnetic form factor
%           is included in the spin-spin correlation function
%           calculation. Possible values:
%               false   No magnetic form factor is applied (default).
%               true    Magnetic form factors are applied, based on the
%                       label string of the magnetic ions, see sw_mff()
%                       function help.
%               cell    Cell type that contains mixed labels and
%                       numbers for every symmetry inequivalent atom in
%                       the unit cell, the numbers are taken as
%                       constants.
%           For example: formfact = {0 'MCr3'}, this won't include
%           correlations on the first atom and using the form factor of
%           Cr3+ ion for the second atom.
% formfactfun
%           Function that calculates the magnetic form factor for given
%           Q value. Default value is @sw_mff(), that uses a tabulated
%           coefficients for the form factor calculation. For
%           anisotropic form factors a user defined function can be
%           written that has the following header:
%               F = @formfactfun(atomLabel,Q)
%           where the parameters are:
%               F   row vector containing the form factor for every
%                   input Q value
%               atomLabel string, label of the selected magnetic atom
%               Q   matrix with dimensions of [3 nQ], where each column
%                   contains a Q vector in Angstrom^-1 units.
% gtensor   If true, the g-tensor will be included in the spin-spin
%           correlation function. Including anisotropic g-tensor or
%           different g-tensor for different ions is only possible
%           here. Including a simple isotropic g-tensor is possible
%           afterwards using the sw_instrument() function.
% useMex    Use mex files to speed up calculations if they are available.
%           default: false
%
% Example:
%
% ...
% horace_on;
% tri = sw_model('triAF',2);
% tri.horace_setpar('mapping',{'J1' 'J2'},'fwhm',0.2);
% d3dobj = d3d(tri.abc,[0 1 0 0],[0,0.01,1],[0 0 1 0],[0,0.01,1],[0 0 0 1],[0,0.1,10]);
% d3dobj = disp2sqw_eval(d3dobj,@cryst.horace,[1 0.5 1.5 0.5 0.01]);
% plot(d3dobj);
%
% This example creates a d3d object, a square in (h,k,0) plane and in
% energy between 0 and 10 meV. Then calculates the inelastice neutron
% scattering intensity of triagular lattice antiferromagnet and plots it
% using sliceomatic.
%
% See also SW, SW.SPINWAVE, SW.MATPARSER, SW.HORACE, SW_READPARAM.
%

inpForm.fname  = {'component' 'norm' 'dE'  'tol' 'optmem' 'hermit' 'notwin'};
inpForm.defval = {'Sperp'     false  0     1e-4  0        true     false};
inpForm.size   = {[1 1]       [1 1]  [1 1] [1 1] [1 1]    [1 1]    [1 1]};
inpForm.soft   = {false       false  false false false    false    false};

inpForm.fname  = [inpForm.fname  {'mapping' 'selector' 'init'}];
inpForm.defval = [inpForm.defval {[]        []         false}];
inpForm.size   = [inpForm.size   {[1 -2]    []         [1 1]}];
inpForm.soft   = [inpForm.soft   {true      true       false}];

inpForm.fname  = [inpForm.fname  {'func'         'formfact' 'formfactfun'}];
inpForm.defval = [inpForm.defval {@obj.matparser false       @sw_mff}];
inpForm.size   = [inpForm.size   {[1 1]          [1 -1]      [1 1]}];
inpForm.soft   = [inpForm.soft   {false          false       false}];

inpForm.fname  = [inpForm.fname  {'gtensor' 'omega_tol' 'useMex'}];
inpForm.defval = [inpForm.defval {false     1e-5        false}];
inpForm.size   = [inpForm.size   {[1 1]     [1 1]       [1 1]}];
inpForm.soft   = [inpForm.soft   {false     false       false}];

warnState = warning('off','sw_readparam:UnreadInput');
param = sw_readparam(inpForm, varargin{:});

% If the field does not exist, create it, else only change parameters
%   specified in the commandline.
if ~isfield(obj.matrix,'horace')
    obj.matrix.horace = param;
else
    for ifl = 1:numel(inpForm.fname)
        idx = find(cellfun(@(c)sum(strcmp(c,inpForm.fname{ifl})),varargin));
        if ~isempty(idx)
            obj.matrix.horace.(inpForm.fname{ifl}) = varargin{idx(end)+1};
        end
    end
end
