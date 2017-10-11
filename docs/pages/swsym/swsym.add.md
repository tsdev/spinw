---
{title: swsym.add, link: swsym.add, summary: saves user defined symmetry operators,
  keywords: sample, sidebar: sw_sidebar, permalink: swsym_add, folder: swsym, mathjax: 'true'}

---

### Syntax

`sym = swsym.add(symstr, {symname})`

### Description

It saves the symmetry generators in symStr into the symmetry.dat file and
returns the line number of the space group in the symmetry.dat file.
 

### Input Arguments

`symStr`
: String, that contains the operators of the space group. If
  not only the generators are given, a possible set of
  generators will be determined and only those will be saved.

`symName`
: Label for the space group.

### See Also

[swsym.generator](swsym_generator) \| [swsym.genreduce](swsym_genreduce)

{% include links.html %}
