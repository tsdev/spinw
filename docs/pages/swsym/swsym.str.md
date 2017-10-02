---
{title: swsym.str, link: swsym.str, summary: generates a string equivalent of symmetry
    operators, keywords: sample, sidebar: sw_sidebar, permalink: swsym_str, folder: swsym,
  mathjax: 'true'}

---

### Syntax

`symstr = swsym.str(symop)`

### Description



### Input Arguments

`symOp`
: Symmetry operator with rotations matrices symOp(:,1:3,:) and
  translation vectors in symOp(:,4,:).

### Output Arguments

strSym    String, contains the symmetry operations.

### See Also

[swsym.add](swsym_add) \| [swsym.generator](swsym_generator)

{% include links.html %}
