---
{title: structfact( ), keywords: sample, summary: calculates magnetic and nuclear structure factor,
  sidebar: sw_sidebar, permalink: '@spinw_structfact.html', folder: '@spinw', mathjax: 'true'}

---
  calculates magnetic and nuclear structure factor
 
  sFact   = STRUCTFACT(obj, kGrid, option1, value1, ...)
 
  sfTable = STRUCTFACT(obj, kGrid, option1, value1, ...)
 
  The calculated structure factors are in barn units. Magnetic structures
  (FM, AFM and HELical) are checked against FullProf. The structure factor
  includes the site occupancy and Debye-Waller factors calculated from
  obj.unit_cell.biso, using the same definition as FullProf.
 
  Input:
 
  obj       Input spinw object, contains crystal and/or magnetic structure.
  kGrid     Defines the reciprocal lattice vectors where the structure
            factor is to be calculated. For commensurate structures these
            are the possible positions of the magnetic Bragg peaks. For
            incommensurate helical/conical structures 3 Bragg peaks
            positions are possible: (k-km,k,k+km) around every reciprocal
            lattice vector. In this case still the integer positions have
            to be given and the code calculates the intensities at all
            three points.
 
  Options:
 
  mode          String, defines the type of calculation:
                    mag     Magnetic structure factor and intensities for
                            unpolarised neutron scattering.
                    nucn    Nuclear structure factor and neutron scattering
                            intensities.
                    nucx    X-ray scattering structure factor and
                            intensities.
  sortq         Sorting the reflections according to increasing momentum
                value if true. Default is false.
 
  gtensor       If true, the g-tensor will be included in the static spin
                correlation function, including anisotropic g-tensor or
                different g-tensor per ion.
 
  formfact      If true, the magnetic form factor is included in the
                spin-spin correlation function calculation. The form factor
                coefficients are stored in obj.unit_cell.ff(1,:,atomIndex).
                Default value is false.
 
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
 
  lambda        Wavelength. If given, the 2theta value for each reflection
                is calculated.
  dmin          Minimum d-value of a reflection, all higher order
                reflections will be removed from the results.
  output        String, defines the type of the output:
                    struct  Results are returned in a struct type variable,
                            default.
                    table   Results are returned in a table type output for
                            easy viewing and exporting.
  tol           Tolerance of the incommensurability of the magnetic
                ordering wavevector. Deviations from integer values of the
                ordering wavevector smaller than the tolerance are considered
                to be commensurate. Default value is 1e-4.
 
  fitmode       Speed up the calculation for fitting mode (omitting
                copying the spinw object to the output). Default is false.
 
  Output:
 
  'sFact' is a structure with the following fields:
  F2            Magnetic structure factor in a matrix with dimensions
                [3 x nHkl].
  Mk            Square of the 3 dimensional magnetic structure factor,
                dimensions are:
                   [nExt(1)*fExt(1) nExt(2)*fExt(2) nExt(3)*fExt(3)],
                where nExt is the size of the extended unit cell.
  hkl           Contains the input Q values, dimensins are [3 nHkl].
  hklA          Same Q values, but in reciproc Angstrom units in the
                lab coordinate system, dimensins are [3 nHkl].
  incomm        Whether the spectra calculated is incommensurate or not.
  formfact      Cell containing the labels of the magnetic ions if form
                factor in included in the spin-spin correlation function.
  {tth}         2theta value of the reflection for the given wavelength,
                only given if a wavelength is provided.
  obj           Copy of the input obj object.
 
  'sfTable' is an optional output in table format for quick viewing and
  saving the output into a text file.
 
  See also SW_QGRID, SW_PLOTSF, SW_INTSF, SPINW.ANNEAL, SPINW.GENMAGSTR.
 
