---
{title: corespec( ), keywords: sample, summary: calculates dispersion and V transformation matrix using linear spin wave theory,
  sidebar: sw_sidebar, permalink: '@spinw_corespec.html', folder: '@spinw', mathjax: 'true'}

---
  calculates dispersion and V transformation matrix using linear spin wave theory
 
  [omega, V] = CORESPEC(obj, hkl, 'option1', value1 ...)
 
  Spin wave dispersion and spin-spin correlation function is calculated at
  the reciprocal space points k. The function can deal with arbitrary
  magnetic structure and magnetic interactions as well as single ion
  anisotropy and magnetic field.
 
  If the magnetic ordering wavevector is non-integer, the dispersion is
  calculated using a coordinate system rotating from cell to cell. In this
  case the Hamiltonian has to fulfill this extra rotational symmetry.
 
  Input:
 
  obj           Input structure, spinw class object.
  hkl           Defines the Q points where the spectra is calculated, in
                reciprocal lattice units, size is [3 nHkl]. Q can be also
                defined by several linear scan in reciprocal space. In this
                case hkl is cell type, where each element of the cell
                defines a point in Q space. Linear scans are assumed
                between consecutive points. Also the number of Q points can
                be specified as a last element, it is 100 by defaults. For
                example: hkl = {[0 0 0] [1 0 0]  50}, defines a scan along
                (h,0,0) from 0 to 1 and 50 Q points are calculated along
                the scan.
 
  Options:
 
  fitmode       Speedup (for fitting mode only), default is false.
  notwin        If true, the spectra of the twins won't be calculated.
                Default is false.
  sortMode      The spin wave modes will be sorted if true. Default is
                true.
  optmem        Parameter to optimise memory usage. The list of hkl values
                will be cut into optmem number of pieces and will be
                calculated piece by piece to decrease memory usage. Default
                of optmem is zero, when the number of slices are determined
                automatically from the available free memory.
  tol           Tolerance of the incommensurability of the magnetic
                ordering wavevector. Deviations from integer values of the
                ordering wavevector smaller than the tolerance are
                considered to be commensurate. Default value is 1e-4.
  omega_tol     Tolerance on the energy difference of degenerate modes when
                diagonalising the quadratic form, default is 1e-5.
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
  nSlice        Parameter to optimise memory usage. The list of hkl values
                will be cut into nSlice pieces and will be calculated piece
                by piece to decrease memory usage. Default is zero, when
                the number of slices are determined automatically from the
                available free memory.
  onlyV         Calculate only the dispersion and the V linear
                transformation matrices. The V matrices transform between
                the original magnon operators and the normal magnon
                operators.
 
  Output:
 
  spectra is a structure, with the following fields:
  omega         Calculated spin wave dispersion, dimensins are
                [nMode nHkl], where nMagExt is the number of magnetic
                atoms in the extended unit cell.
  Sab           Dynamical structure factor, dimensins are
                [3 3 nMode nHkl]. Each (:,:,i,j) submatrix contains the
                9 correlation functions: Sxx, Sxy, Sxz, etc.
  V             Optional output, calculated if 'onlyV' option is true (in
                that case Sab is not calculated). Dimensions are
                [nMode nMode nHkl].
 
 
  nMode is the number of magnetic mode. For commensurate structures it is
  double the number of magnetic atoms in the magnetic cell/supercell. For
  incommensurate structures this number is tripled due to the appearance of
  the (Q+/-km) Fourier components in the correlation functions.
 
  If several twins exist in the sample, omega and Sab are packaged into a
  cell, that contains nTwin number of matrices.
 
  hkl           Contains the input Q values, dimensins are [3 nHkl].
  hklA          Same Q values, but in reciproc Angstrom units in the
                lab coordinate system, dimensins are [3 nHkl].
  incomm        Whether the spectra calculated is incommensurate or not.
  obj           The copy of the input obj.
 
  See also SPINW, SPINW.SPINWAVESYM, SW_NEUTRON, SW_POL, SPINW.POWSPEC, SPINW.OPTMAGSTR.
 
