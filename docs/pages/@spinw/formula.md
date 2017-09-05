---
title: formula( )
keywords: sample
summary: "returns chemical formula, mass, volume, etc."
sidebar: product1_sidebar
permalink: formula.html
folder: @spinw
mathjax: true
---
  returns chemical formula, mass, volume, etc.
 
  formula = FORMULA(obj)
 
  Options:
 
  obj       spinw class object.
 
  Output:
 
  formula struct variable with the following fields:
 
  m             Mass of the unit cell in g/mol unit.
  V             Volume of the unit cell in Angstrom^3 unit.
  rho           Density in g/cm^3 unit.
  chemlabel     List of the different elements.
  chemnum       Number of the listed element names
  chemform      Chemical formula string: series of 'ChemLabel_ChemNum '.
 
  Example:
 
  cryst = spinw('https://goo.gl/do6oTh')
  cryst.formula
 
  The formula of the crystal stored in the test.cif file will be printed
  onto the Command Window.
 
