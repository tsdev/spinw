---
{title: spinw.scga( ), summary: applies the self consistent Gaussian approximation
    at finite temperature, keywords: sample, sidebar: sw_sidebar, permalink: spinw_scga.html,
  folder: spinw, mathjax: 'true'}

---
applies the self consistent Gaussian approximation at finite temperature
 
spectra = SCGA(obj, hkl, 'option1', value1 ...)
 
Input:
 
obj       SpinW object.
hkl       Defines the Q points where the correlations are calculated. It
          is a matrix with dimensions [3,D1,D2,...]. Where the first
          dimesion corresponds to the [h,k,l] index of the Q-point.
 
Options:
 
T         Temperature of the calculation in units given by obj.unit.
plot      If true, the fitting of the integration constant is plotted.
lambda    If given, the integration is avoided.
kbase     Basis vectors that span the Brillouin zone if the system is low
          dimensional. Default value is [] when the dimensionality of the
          system is determined automatically.
nInt      Number of Q points where the Brillouin zone is sampled for the
          integration.
corr      If true, [spinw](spinw.html) correlations are calculated at the given
          momentum vectors. Default is true.
chi       If true, the magnetic susceptibility is calculated.
fitmode       Speedup (for fitting mode only), default is false.
sublat    List of sublattices.
isomode   ...
 
Output:
 
spectra   Structure with fields:
  Sab     Spin-spin correlation function stored in a matrix with
          dimensions of [3,3,1,D1,D2,...].
  Dopt    Optimum value of D.
 

