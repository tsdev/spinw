---
{title: sw_bose, link: sw_bose, summary: coefficient for boson correlation functions,
  keywords: sample, sidebar: sw_sidebar, permalink: sw_bose, folder: swfiles, mathjax: 'true'}

---
  
### Syntax
  
`c = sw_bose(oldt,newt,e)`
  
### Description
  
`c = sw_bose(oldt,newt,e)` calculates the temperature dependent
coefficient for boson correlation functions.
  
### Input Arguments
  
`oldT`
: Original temperature in Kelvin.
  
`newT`
: New temperature in Kelvin.
  
`E`
: Energy in meV, positive is the particle creation side (neutron
  energy loss side in a scattering experiment).
  
### Output Arguments
  
`C`
: Correction coefficients that multiplies the correlation
          function. If any of the input is a vector, `C` will be also a
          vector with the same dimensions.
 

{% include links.html %}
