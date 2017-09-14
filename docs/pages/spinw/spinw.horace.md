---
{title: spinw.horace method, link: spinw.horace, summary: calculates spin wave dispersion/correlation
    functions to be called from Horace, keywords: sample, sidebar: sw_sidebar, permalink: spinw_horace.html,
  folder: spinw, mathjax: 'true'}

---
 
[w, s] = HORACE(obj, qh, qk, ql, 'Option1', Value1, ...)
 
The function produces spin wave dispersion and intensity for Horace
(<a href=http://horace.isis.rl.ac.uk>http://horace.isis.rl.ac.uk</a>).
 
Input:
 
obj           Input spinw object.
qh, qk, ql    Reciprocal lattice components in reciprocal lattice units.
 
Options:
 
component Selects the previously calculated intensity component to be
          convoluted. The possible options are:
              'Sperp' convolutes the magnetic neutron scattering
                      intensity (<Sperp * Sperp> expectation value).
                      Default.
              'Sab'   convolutes the selected components of the spin-spin
                      correlation function. Letter a and b can be 'x',
                      'y' or 'z'. For example: 'Sxx' will convolute the
                      xx component of the correlation function with the
                      dispersion. xyz is the standard coordinate system,
                      see online documentation of SpinW.
          Any linear combination of the above are allowed, for example:
          'Sxx+2*Syy' convolutes the linear combination of the xx
          component of the spin-spin correlation function and the yy
          component.
norm      If true the spin wave intensity is normalized to mbarn/meV/(unit
          cell) units. Default is false.
dE        Energy bin size, for intensity normalization. Use 1 for no
          division by dE in the intensity.
param     Input parameters (can be used also within Tobyfit). Additional
          options ('mat','selector') might be necessary, for details see
          spinw.matparser function. All extra parameters of spinw.horace
          function will be forwarded to the spinw.matparser function before
          calculating the spin wave spectrum (or any user written parser
          function). For user written functions defined with the
          following header:
              func(obj,param)
          the value of the param option will be forwarded. For user
          functions with variable number of arguments, all input options
          of spinw.horace will be forwarded. In this case it is recommended
          to use sw_readparam() function to handle the variable number
          arguments within func().
parfunc   Parser function of the 'param' input. Default is
          @spinw.matparser which can be used directly by Tobyfit. For user
          defined functions the minimum header has to be:
              func(obj,param)
          where obj is an spinw type object, param is the parameter
          values forwarded from spinw.horace directly.
func      User function that will be called after the parameters set on
          the SpinW object. It can be used to optimize magnetic
          structure for the new parameters, etc. The input should be a
          function handle of a function with a header:
              fun(obj)
 
Output:
 
w         Cell that contains the spin wave energies. Every cell elements
          contains a vector of spin wave energies for the corresponding
          input Q vector.
s         Cell that contains the calculated element of the spin-spin
          correlation function. Every cell element contains a vector of
          intensities in the same order as the spin wave energies in w.
 
Example:
 
...
horace_on;
d3dobj = d3d(cryst.abc,[0 1 0 0],[0,0.01,1],[0 0 1 0],[0,0.01,1],[0 0 0 1],[0,0.1,10]);
d3dobj = disp2sqw_eval(d3dobj,@cryst.horace,{'component','Sperp'},0.1);
plot(d3dobj);
 
This example creates a d3d object, a square in (h,k,0) plane and in
energy between 0 and 10 meV. Then calculates the inelastice neutron
scattering intensity of the spin wave model stored in cryst and plots it
using sliceomatic.
 
See also SPINW, SPINW.SPINWAVE, SPINW.MATPARSER, SW_READPARAM.
 

