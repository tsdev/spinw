---
{title: sw_neutron, link: sw_neutron, summary: calculates neutron scattering cross
    section, keywords: sample, sidebar: sw_sidebar, permalink: sw_neutron, folder: swfiles,
  mathjax: 'true'}

---
  
### Syntax
  
`spectra = sw_neutron(spectra,Name,Value)`
  
### Description
  
`spectra = sw_neutron(spectra,Name,Value)` calculates the neutron
scattering cross section for polarised and unpolarised neutrons. The
function reads the calculated spin-spin correlation function
$$\mathcal{S}^{\alpha\beta}(\mathbf{Q},\omega)$$ and calculates the neutron
scattering cross section for unpolarized neutrons using the formula:
 
$$S_\perp(Q,\omega)=\sum_{\alpha\beta}(1-\hat{q}^\alpha\hat{q}^\beta)\cdot S^{\alpha\beta}(Q,\omega)$$
   
It also calculates spin-spin correlation function in the Blume-Maleev
coordinate system and the complete polarised neutron scattering cross
section.
  
{% include note.html content=" The Blume-Maleev coordinate system is a cartesian coordinate
system with $$x_{BM}$$, $$y_{BM}$$ and $$z_{BM}$$ basis vectors defined as:
<br> $$x_{BM}$$    parallel to the momentum transfer $$Q$$,
<br> $$y_{BM}$$    perpendicular to $$x_{BM}$$ in the scattering plane,
<br> $$z_{BM}$$    perpendicular to the scattering plane.
" %}
 
### Input Arguments
  
`spectra`
: Input structure, contains spin-spin correlation functions. Supported
  inputs are produced by [spinw.spinwave](spinw_spinwave), [spinw.powspec](spinw_powspec) and
  [spinw.scga].
  
### Name-Value Pair Arguments
  
`'n'`
: Normal vector to the scattering plane, in real space ($$xyz$$
  coordinate system), stored in a row vector with 3 elements. Default
  value is `[0 0 1]`.
 
`'uv'`
: Cell, that contains two vectors defining the scattering 
  plane in rlu. If given overwrites the `n` parameter value. For example:
  `{[1 0 0] [0 1 0]}` stands for the $$(h,k,0)$$ scattering plane.
  
`'pol'`
: If `true` the cross sections in the Blume-Maleev
  coordinate system will be also calculated (`inP`, `Pab` and `Mab`
  fields of the output `spectra`). Default value is `false`.
  
### Output Arguments
  
`spectra`
: Same as the input `spectra` plus the following additional fields:
  * `param`   Input parameters.
  * `Sperp`   $$S_\perp(i_{mode},\mathbf{Q})$$ unpolarised neutron 
              scattering cross section, stored in a matrix with
              dimensions of $$[n_{mode}\times n_{hkl}]$$.
  * `intP`    $$I_P(P_i,i_{mode},\mathbf{Q})$$ polarised neutron scattering 
              cross section, when only the incident neutron polarization
              is analyzed. It is stored in a matrix with dimensions of
              $$[3\times n_{mode}\times n_{hkl}]$$.
  * `Pab`     $$I_{Pab}(P_f,P_i,i_{mode},\mathbf{Q})$$ complete polarised 
              neutron scattering cross section, when the polarisation of
              both the incident ($$P_i$$) and the scattered ($$P_f$$)
              neutrons are analyzed. Stored in a matrix with dimensions
              of $$[3\times 3\times n_{mode}\times n_{hkl}]$$.
  * `Mab`     $$M_{ab}(P_f,P_i,i_{mode},\mathbf{Q})$$ components of the 
              spin-spin correlation function in the blume-Maleev
              coordinate system, stored in a matrix with dimensions of
              $$[3\times \times3 n_{mode}\times n_{hkl}]$$.
 
If several domains exist in the sample, `Sperp`, `intP`, `Pab` and `Mab`
will be packaged into a cell, that contains $$n_{twin}$$ number of
matrices.
 
The meaning of the indices above:
* $$P_i$$: index of incident polarisation ($$1=x_{BM}$$, $$2=y_{BM}$$ or $$3=z_{BM}$$),
* $$P_f$$: index of final polarisation ($$1=x_{BM}$$, $$2=y_{BM}$$ or $$3=z_{BM}$$),
* $$i_{mode}$$: index of spin wave mode,
* $$\mathbf{Q}$$: index of momentum transfer.
 
  
### See Also
  
[sw_egrid](sw_egrid) \| [spinw](spinw) \| [spinw.spinwave](spinw_spinwave)
 
*[rlu]: reciprocal lattice units
 

{% include links.html %}
