---
{title: sw_model, link: sw_model, summary: creates predefined spin models, keywords: sample,
  sidebar: sw_sidebar, permalink: sw_model, folder: swfiles, mathjax: 'true'}

---
  
### Syntax
  
`obj = sw_model(model, param)`
  
### Description
  
`obj = sw_model(model, param)` generates spin models, such as triangular
lattice antiferromagnet, square lattice, etc. It also generates the
magnetic ground state. For each lattice arbitrary number of further
neighbor bonds can be defined using a vector of exchange values.
  
### Input Arguments
  
`model`
: String, name of the model, one of the following:
  * `'triAF'`     Triangular lattice Heisenberg antiferromagnet
                  in the $$ab$$ plane ($$a=b=3$$ Å), with γ =
                  120° angle and optimised magnetic structure.
  * `'squareAF'`  Square lattice antiferromagnet.
  * `'chain'`     Chain with further neighbor interactions.
  
`param`
: Input parameters of the model, row vector which gives the values of the
  Heisenberg exchange for first, second, thirs etc. neighbor bonds stored
  in `p(1)`, `p(2)`, `p(3)`, etc. respectively.
  
### Output Arguments
  
`obj`
: [spinw](spinw) class object with the selected model.
  
### See Also
  
[spinw](spinw)
 

{% include links.html %}
