---
{title: sw_econtract( ), summary: 'converts (Q,omega) values to Qm values for diffraction
    instrument', keywords: sample, sidebar: sw_sidebar, permalink: swfiles_sw_econtract.html,
  folder: swfiles, mathjax: 'true'}

---
converts (Q,omega) values to Qm values for diffraction instrument
 
Qm = SW_FLATCONE(Q,'Option1', Value,...) 
 
Input:
 
Q         Input values in reciprocal space in the scattering plane in
          Angstrom^-1 units, dimensions are [2 nQ].
 
Options:
 
omega     Energy transfer value in meV, default is zero.
lambda    Wavelength of the incident beam in Angstrom.
ki        Momentum of the incidend beam in Angstrom^-1.
sense     Scattering sense:
              1       detectors are on the right hand side from the
                      incident beam direction, default.
             -1       detectors are on the left hand side from the
                      incident beam direction.
 
See also SW_CONVERTER.
 
