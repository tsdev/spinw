---
title: str( )
keywords: sample
summary: "generates a string equivalent of symmetry operators"
sidebar: product1_sidebar
permalink: str.html
folder: +swsym
mathjax: true
---
  generates a string equivalent of symmetry operators
 
  symStr = SWSYM.STR(symOp)
 
  Input:
 
  symOp     Symmetry operator with rotations matrices symOp(:,1:3,:) and
            translation vectors in symOp(:,4,:).
 
  Output:
 
  strSym    String, contains the symmetry operations.
 
  See also SWSYM.ADD, SWSYM.GENERATOR.
 
