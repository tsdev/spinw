---
{title: spinw.symbolic method, link: spinw.symbolic, summary: switches between symbolic/numeric
    mode, keywords: sample, sidebar: sw_sidebar, permalink: spinw_symbolic, folder: spinw,
  mathjax: 'true'}

---
  
### Syntax
  
`symb = symbolic(obj)`
 
`symbolic(obj, symb)`
  
### Description
  
`symb = symbolic(obj)` returns `true` if symbolic calculation mode is on,
`false` for numeric mode.
   
`symbolic(obj, symb)` sets whether the calculations are in
symbolic/numeric (`true`/`false`) mode. Switching to symbolic mode, the
spin values, matrix elements, magnetic field, magnetic structure and
physical units are converted into symbolic variables. If this is not
desired, start with a symbolic mode from the beggining and have full
control over the values of the above mentioned variables.
  
### See Also
  
[spinw](spinw) \| [spinw.spinwavesym](spinw_spinwavesym)
 

{% include links.html %}
