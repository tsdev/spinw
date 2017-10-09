---
{title: spinw.table method, link: spinw.table, summary: outputs easy to read tables
    of internal data, keywords: sample, sidebar: sw_sidebar, permalink: spinw_table,
  folder: spinw, mathjax: 'true'}

---
  
### Syntax
  
`T = table(obj,type,{index},{showval})`
  
### Description
  
`T = table(obj,type,{index},{showval})` returns a table that shows in an
easy to read/export format different internal data, such as magnetic atom
list, bond list, magnetic structure, etc.
   
For the matrix labels in the list of bonds, the '>>' sign means that the
matrix value is determined using the bond symmetry operators.
   
{% include note.html content=" The `table` data type is only supported in Matlab R2013b or newer.
When running older versions of Matlab, `spinw.table` returns a struct." %}
  
### Input Arguments
  
`obj`
: [spinw](spinw) object.
  
`type`
: String, determines the type of data to show, possible values are:
  * `'matom'`     properties of magnetic atoms in the unit cell,
  * `'matrix'`    list of matrices,
  * `'ion'`       single ion term in the Hamiltonian,
  * `'bond'`      properties of selected bonds,
  * `'mag'`       magnetic structure.
  
`index`
: Indexing into the type of data to show, its meaning depends on the
  `type` parameter. For `'bond'` indexes the bonds (1 for first
  neighbors, etc.), if empty all bonds will be shown. For `'mag'` it
  indexes the propagation vectors, the magnetization of the selected
  propagation vector will be shown. Default value is 1, if empty vector `[]` is given, all
  bonds/propagation vector will be shown.
  
`showVal`
: If `true`, also the values of the single ion and exchange matrices
  will be shown. The values shown  are the symmetry transformed exchange
  values after the symmetry operations (if there is any). Default value
  is `false`.
  
### Output Arguments
  
`T`
: Matlab `table` type object.
 

{% include links.html %}
