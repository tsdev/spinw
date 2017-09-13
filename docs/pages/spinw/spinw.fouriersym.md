---
{title: spinw.fouriersym( ), summary: calculates the Fourier transformation of a symbolic
    Hamiltonian, keywords: sample, sidebar: sw_sidebar, permalink: spinw_fouriersym.html,
  folder: spinw, mathjax: 'true'}

---
calculates the Fourier transformation of a symbolic Hamiltonian
 
res = FOURIER(obj, 'option1', value1 ...)
 
Input:
 
obj           Input structure, spinw class object.
 
Options:
 
hkl           Symbolic definition of q vector. Default is the general Q
              point:
                  hkl = [sym('h') sym('k') sym('l')]
 
 
 
See also SPINW.FOURIER.
 
