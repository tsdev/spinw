---
{title: sw_magdomain, link: sw_magdomain, summary: calculates the spin-spin correlation
    function for magnetic domains, keywords: sample, sidebar: sw_sidebar, permalink: sw_magdomain,
  folder: swfiles, mathjax: 'true'}

---
  
### Syntax
  
`spectra = sw_magdomain(spectra,Name,Value)`
  
### Description
  
`spectra = sw_magdomain(spectra,Name,Value)` calculates spin-spin
correlation function averaged over magnetic domains that are related by a
point group operation. Several domains with different volume ratios can
be defined. The spin-spin correlation function will be rotated and summed
according to the domains. The rotations of the magnetic domains are
defined in the $$xyz$$ coordinate system, same as the coordinate system for
the spin-spin correlation function. The function only rotates the
$$\mathcal{S}^{\alpha\beta}$$ components of the spin, but not the momentum
$$Q$$, thus it cannot be used to simulate magnetic domains with different
propagation vector.
  
### Examples
  
The above example calculates the spectrum for magnetic domains that are
related by a 90 ° rotation around the $$z$$-axis (perpendicular to the
$$ab$$ plane). All domains have equal volume.
  
```matlab
spec = cryst.spinwave({[0 0 0] [1 0 0]})
spec = sw_magdomain(spec,'axis',[0 0 1],'angled',[0 90 180 270]);
```
  
### Input Arguments
  
`spectra`
: Calculated spin wave spectrum.
  
### Name-Value Pair Arguments
  
`'axis'`
: Defines axis of rotation to generate domains in the $$xyz$$
  coordinate system, row vector with 3 elements.
  
`'angle'`
: Defines the angle of rotation to generate domains in radian
  units, multiple domains can be defined if angle is a
  row vector with $$n_{dom}$$ number of elements.
  
`'angled'`
: Same as the `angle` parameter, just in ° units.
  
`'rotC'`
: Rotation matrices, that define crystallographic domains, alternative
  input instead of `angle` and `axis`, matrix with dimensions of
  $$[3\times 3\times n_{dom}]$$.
  
`'vol'`
: Volume fractions of the domains in a row vector with $$n_{dom}$$ number of
  elements. Default value is `ones(1,nDom)`.
  
### Output Arguments
  
`spectra`
: Spectrum (Struct) with the following additional fields:
  * `Sab`     The multi domain spectrum will be stored here.
  * `Sabraw`  The original single domain spectrum is kept here, so that a
              consecutive run of `sw_magdomain` will use the original single
              domain spectrum, without the need of recalculating the full
              spectrum.
  * `domVol`  Volume of each domains in a row vector with $$n_{dom}$$
              number of elements.
  * `domRotC` Rotation matrices for each domain, with dimensions of
              $$[3\times 3\times n_{dom}]$$.
  
### See Also
  
[spinw.spinwave](spinw_spinwave) \| [spinw.addtwin](spinw_addtwin) \| [spinw.twinq](spinw_twinq)
 

{% include links.html %}
