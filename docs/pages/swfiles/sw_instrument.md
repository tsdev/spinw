---
{title: sw_instrument( ), summary: includes instrumental factors into the calculated
    spectrum, keywords: sample, sidebar: sw_sidebar, permalink: sw_instrument.html,
  folder: swfiles, mathjax: 'true'}

---
includes instrumental factors into the calculated spectrum
 
spectra = SW_INSTRUMENT(spectra, 'option1','value1',...)
 
It includes instrumental factors (resolution, energy transfer range,
etc.) to the convoluted spectrum.
 
Options:
 
dE            Defines the FWHM energy resolution of the instrument. It
              can be a string, single number, vector of function hangle:
                String    File name, that contains the FWHM energy
                          resolution values as a function of energy
                          transfer. The file has to contain two columns,
                          first is the energy values, the second is the
                          FWHM resolution at the given energy transfer
                          value, see [sw_res()](sw_res.html) function for details.
                Number    Constant FWHM energy resolution as a function
                          of energy transfer.
                Matrix    Dimensions of Nx2, first column contains the
                          energy transfer values, second column contains
                          the FWHM resolution values. These values will
                          be fitted using a polynomial with a fixed
                          degree, see [sw_res()](sw_res.html) function for details.
                Function  Function handle of a resolution function
                          with the following header:
                              E_FWHM = res_fun(E)
                          where E_FWHM is the FWHM energy resolution and
                          E is the energy transfer value.
func          Shape of the energy resolution function, for details see
              the help of [sw_resconv](sw_resconv.html).
polDeg        Degree of the fitted polynomial to the instrumental
              resolution data. Default is 5.
dQ            Momentum transfer resolution of the instrument, FWHM is
              given in A-1 units, default is 0.
ThetaMin      Minimum scattering angle in degree, default is 0.
plot          If the resolution is read from file and plot option is
              true, tre resolution will be plotted, default is true.
 
Fixed incident neutron energy:
ki            Momentum of the incident neutrons in A^-1 units.
Ei            Energy of the incident neutrons in meV.
 
Fixed final neutron energy:
kf            Final momentum of the neutrons in A^-1 units.
Ef            Final neutron energy in meV.
 
norm          If true, the data is normalized to mbarn units. Default is
              false. If no g-tensor is included in the spin wave
              calculation, g-tensor = 2 is assumed here.
useRaw        If false, the already modified spectra.swConv field is
              modified further instead of the original powder spectrum
              stored in spectra.swRaw. Default is true.
 
Output:
 
spectra       Struct variable, same as input with following additional
              fields:
 
norm          True, if the spectrum is normalised to mbarn units.
ki            Incident neutron wave vector as given in the input.
dE            Energy resolution polynomial as given in the input.
dQ            FWHM of the momentum resolution.
swRaw         Original simulated data, withouth the application of the
              instrumental factors.
 
 
See also POLYFIT, POLYVAL, SW_RES, SW_RESCONV.
 

