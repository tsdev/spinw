---
{title: swsym.genreduce, link: swsym.genreduce, summary: reduces symmetry operators
    to the generators, keywords: sample, sidebar: sw_sidebar, permalink: swsym_genreduce.html,
  folder: swsym, mathjax: 'true'}

---
  
### Syntax
  
`[symOpG, isGen] = swsym.genreduce(symOp)`
  
### Description
  
`[symOpG, isGen] = swsym.genreduce(symOp)` takes the list of symmetry
operators in `symOp` and determines a minimum subset of operators that
can generate the given list.
  
### Input Arguments
  
`symOp`
: Matrix that contains both the rotation and translation matrices
  having dimensions of $$[3\times 4\times n_{sym}]$$, where the
  `symMat(:,4,:)` stores the translation vectors, while the
  `symMat(:,1:3,:)` stores the $$3\times 3$$ rotation matrices.
  
### Output Arguments
  
`symOpG`
: A set of operators, that can generate all the operators of the input.
 
`isGen`
: Vector, that gives whether a given input operator is part of
  the generators, dimensions are $$[1\times n_{sym}]$$.
  
### See Also
  
[swsym.add](swsym_add.html) \| [swsym.generator](swsym_generator.html) \| [swsym.operator](swsym_operator.html)
 

{% include links.html %}
