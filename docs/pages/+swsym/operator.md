---
title: operator( )
keywords: sample
summary: "calculates all symmetry operators or general positions for a space group"
sidebar: product1_sidebar
permalink: operator.html
folder: +swsym
mathjax: true
---
  calculates all symmetry operators or general positions for a space group
 
  [symOp, symInfo] = SWSYM.OPERATOR(sym, fid)
 
  Input:
 
  sym           Line index in the symmetry.dat file or string of the
                symmetry operators or matrix of symmetry generators with
                dimensions of [3 4 nOp].
                For example:
                    sym = 'P b n m';
  fid           For printing the symmetry operators:
                    0   no printed output (Default)
                    1   standard output (Command Line)
                    fid text file opened before with the fid = fopen(path)
 
  Output:
 
  symOp         The rotational part of the symmetry operators, dimensions
                are [3 3 nSym].
  symInfo       Structure containing additional information about the space
                group with the following fields:
    name            Name of the space group in string. If function called
                    with no input, name stores the name of all spase groups
                    from symmetry.dat in a cell.
    str             The string of the symmetry operations.
    num             The index of the symmetry in the symmetry.dat file.
 
 
  See also SPINW, SWSYM.POSITION, SPINW.ATOM, SPINW.MATOM, SWSYM.GENERATOR,
  SWSYM.POINT.
 
