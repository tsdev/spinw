---
{title: sw_atomdata, link: sw_atomdata, summary: returns information on chemical elements,
  keywords: sample, sidebar: sw_sidebar, permalink: sw_atomdata, folder: swfiles,
  mathjax: 'true'}

---
  
### Syntax
  
`data = sw_atomdata(atomsymb)`
  
`data = sw_atomdata(atomsymb,datatype)`
 
### Description
  
`data = sw_atomdata(atomsymb)` returns information on chemical elements
(RGB color code, mass, long name) in a struct. The element is identified
by its short name, such as 'O' for oxygen. If the given atom name does
not exists, the function returns the data for `'Unobtanium'`.
  
`data = sw_atomdata(atomsymb,datatype)` returns only the requested type
of data.
 
### Examples
 
The radius of the hydrogen atom:
```matlab
r_H = sw_atomdata('H','radius')
```
*Output*
```
r_H =
    0.3700
```
 
  
### Input Arguments
  
`atomSymb`
: String of the name of the atom, for example 'He' or the atomic
  number Z. If the string contains whitespace character, the
  second word will be used to identify the atom.
  
`dataType`
: Type of information requested, following strings are accepted:
  * `'name'`        Atomic symbol.
  * `'radius'`      Atomic radius.
  * `'color'`       Color of the atom from the [CPK color scheme](https://en.wikipedia.org/wiki/CPK_coloring).
  * `'mass'`        Average mass of the element.
  * `'longname'`    Name of the element.
  * `'Z'`           Atomic index.
  * `'all'`         All atomic data returned in a struct.
  
### See Also
  
[sw_mff](sw_mff) \| [sw_cff](sw_cff)
 

{% include links.html %}
