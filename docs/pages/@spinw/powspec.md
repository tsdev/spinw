---
title: powspec( )
keywords: sample
summary: "calculates powder averaged spin wave spectra"
sidebar: product1_sidebar
permalink: powspec.html
folder: @spinw
mathjax: true
---
  calculates powder averaged spin wave spectra
 
  spectra = POWSPEC(obj, hklA, 'Option1', Value1, ...)
 
  The function calculates powder averaged spectrum by doing a 3D average in
  momentum space. This method is not efficient for low dimensional (2D, 1D)
  structures. To speed up the calculation with mex files use the
  swpref.setpref('usemex',true) option. The function can do powder average
  on arbitrary spectral function, but it is currently tested with two
  functions:
        spinw.spinwave  Powder average spin wave spectrum.
        spinw.scga      Powder averaged diffuse scattering spectrum.
  The type of spectral function is determined by the specfun option.
 
  Input:
 
  obj       spinw class object.
  hklA      Vector containing the Q values in inverse Angstrom where powder
            spectra will be calculated, dimensions are [1 nQ].
 
  Options:
 
  specfun   Function handle of the spectrum calculation function. Default
            is @spinwave.
  nRand     Number of random orientations per Q value, default is 100.
  Evect     Vector, defines the center/edge of the energy bins of the
            calculated output, dimensions are is [1 nE]. The energy units
            are defined by the unit.kB property of the spinw object. Default
            value is an edge bin: linspace(0,1.1,101).
  binType   String, determines the type of bin give, possible options:
                'cbin'    Center bin, the center of each energy bin is given.
                'ebin'    Edge bin, the edges of each bin is given.
            Default is 'ebin'.
  T         Temperature to calculate the Bose factor in units
            depending on the Boltzmann constant. Default is taken from
            obj.single_ion.T value.
  title     Gives a title string to the simulation that is saved in the
            output.
  extrap    If true, arbitrary additional parameters are passed over to
            the spectrum calculation function.
  fibo      If true, instead of random sampling of the unit sphere the
            Fibonacci numerical integration is implemented as described in:
            J. Phys. A: Math. Gen. 37 (2004) 11591
            The number of points on the sphere is given by the largest
            Fibonacci number below nRand. Default is false.
  imagChk   Checks that the imaginary part of the spin wave dispersion is
            smaller than the energy bin size. Default is true.
  component See the help of sw_egrid() function for description.
 
  The function accepts all options of spinw.spinwave() with the most
  important options are:
 
  formfact      If true, the magnetic form factor is included in the
                spin-spin correlation function calculation. Default value
                is false.
  formfactfun   Function that calculates the magnetic form factor for given
                Q value. Default value is @sw_mff(), that uses a tabulated
                coefficients for the form factor calculation. For
                anisotropic form factors a user defined function can be
                written that has the following header:
                    F = @formfactfun(atomLabel,Q)
                where the parameters are:
                    F   row vector containing the form factor for every
                        input Q value
                    atomLabel string, label of the selected magnetic atom
                    Q   matrix with dimensions of [3 nQ], where each column
                        contains a Q vector in Angstrom^-1 units.
  gtensor       If true, the g-tensor will be included in the spin-spin
                correlation function. Including anisotropic g-tensor or
                different g-tensor for different ions is only possible
                here.
  hermit        Method for matrix diagonalization:
                    true      J.H.P. Colpa, Physica 93A (1978) 327,
                    false     R.M. White, PR 139 (1965) A450.
                Colpa: the grand dynamical matrix is converted into another
                       Hermitian matrix, that will give the real
                       eigenvalues.
                White: the non-Hermitian g*H matrix will be diagonalised,
                       that is not the elegant method.
                Advise:
                Always use Colpa's method, except when small imaginary
                eigenvalues are expected. In this case only White's method
                work. The solution in this case is wrong, however by
                examining the eigenvalues it can give a hint where the
                problem is.
                Default is true.
 
  The function accepts some options of spinw.scga() with the most important
  options are:
 
  nInt      Number of Q points where the Brillouin zone is sampled for the
            integration.
 
  Output:
 
  'spectra' is a struct type variable with the following fields:
  swConv    The spectra convoluted with the dispersion. The center
            of the energy bins are stored in spectra.Evect. Dimensions are
            [nE nQ].
  hklA      Same Q values as the input hklA [1 nQ]. Evect
            Contains the input energy transfer values, dimensions are
            [1 nE].
  param     Contains all the input parameters.
  obj       The copy of the input obj object.
  Evect     Energy grid converted to edge bins.
 
  Example:
 
  tri = sw_model('triAF',1);
  E = linspace(0,4,100);
  Q = linspace(0,4,300);
  triSpec = tri.powspec(Q,'Evect',E,'nRand',1e3);
  sw_plotspec(triSpec);
 
  The example calculates the powder spectrum of the triangular lattice
  antiferromagnet (S=1, J=1) between Q = 0 and 3 A^-1 (the lattice
  parameter is 3 Angstrom).
 
  See also SPINW, SPINW.SPINWAVE, SPINW.OPTMAGSTR, SW_EGRID.
 
