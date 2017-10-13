---
{title: sw_cartesian, link: sw_cartesian, summary: creates a right handed Cartesian
    coordinate system, keywords: sample, sidebar: sw_sidebar, permalink: sw_cartesian,
  folder: swfiles, mathjax: 'true'}

---
  
### Syntax
  
`[vy, vz, vx] = sw_cartesian(n)`
 
`V = sw_cartesian(n)`
  
### Description
  
`[vy, vz, vx] = sw_cartesian(n)` creates an $$(x,y,z)$$ right handed
Cartesian coordinate system with $$v_x$$, $$v_y$$ and $$v_z$$ defining the
basis vectors. 
 
`V = sw_cartesian(n)` the generated basis vectors are stored in the `V`
matrix: `V = [vx vy vz]` as column vectors.
  
### Input Arguments
  
`n`
: Either a 3 element row/column vector or a $$[3\times 3]$$ matrix with
  columns defining 3 vectors.
  
### Output Arguments
  
`vy,vz,vx`
: Vectors defining the right handed coordinate system. They are
          either column of row vectors depending on the shape of the
          input `n`.
 

{% include links.html %}
