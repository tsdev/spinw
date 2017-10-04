---
{title: spinw.field method, link: spinw.field, summary: get/set magnetic field value,
  keywords: sample, sidebar: sw_sidebar, permalink: spinw_field, folder: spinw, mathjax: 'true'}

---
  
### Syntax
  
`field(obj,B)`
`B = field(obj)`
  
### Description
  
`field(obj,B)` sets the magnetic field stored in `obj.single_ion.field`
to `B`, where `B` is a $$[1\times 3]$$ vector.
   
`B = field(obj)` returns the current value of the magnetic field value
stored in `obj`.
   
### See Also
  
[spinw](spinw) \| [spinw.temperature](spinw_temperature) \| [spinw.single_ion](spinw_single_ion)
 

{% include links.html %}
