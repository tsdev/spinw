---
{title: sw_neutron( ), keywords: sample, summary: calculates neutron scattering intensity for spin wave spectrum,
  sidebar: sw_sidebar, permalink: swfiles_sw_neutron.html, folder: swfiles, mathjax: 'true'}

---
  calculates neutron scattering intensity for spin wave spectrum
 
  spectra = SW_NEUTRON(spectra, 'Option1', Value1' ...)
 
  It calculates neutron scattering intensity for polarised and unpolarised
  neutron scattering on spin waves.
 
  Input:
 
  spectra   Input structure that contains the spin-spin correlation
            function.
 
  Options:
 
  uv        Cell, that contains two vectors, that define the scattering 
            plane in r.l.u. For example: {[1 0 0] [0 1 0]} for the hk
            plane.
  n         Normal vector to the scattering plane, in real space (xyz
            coordinate system), dimensions are [1 3]. Default is [0 0 1].
  pol       Whether to calculate cross sections in the Blume-Maleev
            coordinate system (inP, Pab and Mab fields of spectra). Default
            is false.
 
  Output:
 
  the spectra output has the following additional fields:
  param     Input parameters.
 
  Sperp     Sperp(mode,Q) unpolarised neutron scattering cross section,
            dimensions are [nMode nHkl].
  intP      intP(Pi,mode,Q) polarised scattering cross section, dimensions
            are [3 nMode nHkl].
  Pab       Pab(Pf,Pi,mode,Q) complete polarised scattering cross section,
            dimensions are [3 3 nMode nHkl].
  Mab       Mab(Pf,Pi,mode,Q) components of the correlation function,
            dimensions are [3 3 nMode nHkl].
 
  If several domains exist in the sample, Sperp, intP, Pab and Mab are
  packaged into a cell, that contains nTwin number of matrices.
 
  The meaning of the indices:
            Pi      index of incident polarisation (1=xBM, 2=yBM or 3=zBM),
            Pf      index of final polarisation (1=xBM, 2=yBM or 3=zBM),
            mode    index of spin wave mode,
            Q       index of momentum transfer.
 
  The polarisation components Pi and Pf defines the Blume-Maleev coordinate
  system with (xBM, yBM and zBM) basis vectors as follows:
            xBM     parallel to the momentum transfer Q,
            yBM     perpendicular to Px in the scattering plane,
            zBM     perpendicular to the scattering plane.
 
  See also SPINW, SPINW.SPINWAVE.
 
