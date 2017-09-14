---
{title: swsym.str( ), summary: generates a string equivalent of symmetry operators,
  keywords: sample, sidebar: sw_sidebar, permalink: swsym_str.html, folder: swsym,
  mathjax: 'true'}

---
 
symStr = SWSYM.STR(symOp)
 
Input:
 
symOp     Symmetry operator with rotations matrices symOp(:,1:3,:) and
          translation vectors in symOp(:,4,:).
 
Output:
 
strSym    String, contains the symmetry operations.
 
See also SWSYM.ADD, SWSYM.GENERATOR.
 

