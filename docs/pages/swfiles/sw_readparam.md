---
title: sw_readparam( )
keywords: sample
summary: "parse input arguments (option, value pairs)"
sidebar: product1_sidebar
permalink: sw_readparam.html
folder: swfiles
mathjax: true
---
  parse input arguments (option, value pairs)
 
  input = SW_READPARAM(format, raw)
 
  It reads in parameters from input structure. Lower and upper case
  insensitive, the output structure has the names stored in format.fname.
  Instead of a struct type input, also a list of parmeters can be given in
  a parameter name, value pairs. Where the parameter name is a string.
 
  Input:
 
  format is struct type with the following fields:
  fname     Field names, strings in cell, dimensions are [nParm 1].
  size      field size, if negative means index, field sizes with same
            negative index have to be the same size.
  defval    Optional, default value if missing.
  soft      Optional, if exist and equal to 1, in case of bad input
            value, defval is used without error message.
 
