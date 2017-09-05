---
title: sw_atomdata( )
keywords: sample
summary: "returns information on elements stored in the atom.dat file"
sidebar: product1_sidebar
permalink: sw_atomdata.html
folder: swfiles
mathjax: true
---
  returns information on elements stored in the atom.dat file
 
  [data atomLabel] = SW_ATOMDATA(atomSymb, {dataType})
 
  Input:
 
  atomSymb  String of the name of the atom, for example 'He' or the atomic
            number Z. If the string contains whitespace character, the
            second word will be used to identify the atom.
  dataType  Type of information requested:
                name        Atomic symbol.
                radius      Atomic radius.
                color       Color of the atom from the CPK color scheme.
                mass        Average mass of the element.
                longname    Name of the element.
                Z           tomic index.
                all         All atomic data returned in a struct.
 
  Example:
  sw_atomdata('H','radius') = 0.37
  If the atom label does not exists, the function returns radius = 1,
  color = [255 167 0].
 
  optional second output is 'atomLabel' that contains the name of the atom
  clean.
 
  See also SW_MFF, SW_CFF.
 
