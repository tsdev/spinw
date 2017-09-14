---
{title: swsym.genreduce( ), summary: reduces the list of symmetry operators to the
    generators, keywords: sample, sidebar: sw_sidebar, permalink: swsym_genreduce.html,
  folder: swsym, mathjax: 'true'}

---
 
[symOpG, isGen] = SWSYM.GENREDUCE(symOp)
 
Input:
 
symOp     Matrix that contains both the rotation and translation matrices
          having dimensions of [3 4 nSym], where the symMat(:,4,:) stores
          the translation vectors, while the symMat(:,1:3,:) stores the
          rotation operators.
 
Output:
 
symOpG    A set of operators, that can generate all the operators of the
          input.
isGen     Vector, that gives whether a given input operator is part of
          the generators, dimensions are [1 nSym].
 
See also SWSYM.ADD, SWSYM.GENERATOR, SWSYM.OPERATOR.
 

