---
{title: spinw.fouriersym method, link: spinw.fouriersym, summary: calculates the Fourier
    transformation of the symbolic Hamiltonian, keywords: sample, sidebar: sw_sidebar,
  permalink: spinw_fouriersym, folder: spinw, mathjax: 'true'}

---
  
### Syntax
  
`res = fouriersym(obj,Name,Value)`
  
### Description
  
`res = fouriersym(obj,Name,Value)` solves the symbolic Fourier transform
problem.
  
### Input Arguments
  
`obj`
: [spinw](spinw) object.
  
### Name-Value Pair Arguments
  
`'hkl'`
: Symbolic definition of positions in momentum space. Default value is
  the general $$Q$$ point:
  ```matlab
  hkl = [sym('h') sym('k') sym('l')]
  ```
  
### See Also
  
[spinw.fourier](spinw_fourier)
 

{% include links.html %}
