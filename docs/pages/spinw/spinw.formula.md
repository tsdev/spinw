---
{title: spinw.formula method, link: spinw.formula, summary: returns basic physical
    properties, keywords: sample, sidebar: sw_sidebar, permalink: spinw_formula, folder: spinw,
  mathjax: 'true'}

---
  
### Syntax
  
`formula = formula(obj)`
  
### Description
  
`result = formula(obj)` returns chemical mass, density, cellvolume etc.
of `obj`.
  
### Examples
  
The formula of the crystal stored in the
[https://goo.gl/do6oTh](https://goo.gl/do6oTh) linked file will be
printed onto the Command Window.
 
```matlab
cryst = spinw('https://goo.gl/do6oTh')
cryst.formula
```
*Output*
```
     Chemical formula:  Cr1Li1O2
     Formula mass:        90.936 g/mol
     Formula in cell:          3 units
     Cell volume:        105.178 Angstrom^3
     Density:              4.307 g/cm^3
```
 
  
### Name-Value Pair Arguments
  
`'obj'`
: [spinw](spinw) object.
  
### Output Arguments
  
`formula` struct variable with the following fields:
* `m`         Mass of the unit cell in g/mol units.
* `V`         Volume of the unit cell in Å$$^3$$ units.
* `rho`       Density in g/cm$$^3$$ unit.
* `chemlabel` List of the different elements.
* `chemnum`   Number of the listed element names
* `chemform`  Chemical formula string.
 

{% include links.html %}
