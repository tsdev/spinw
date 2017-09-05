---
title: sw_cff( )
keywords: sample
summary: "returns the atomic charge form factor values for X-ray scattering"
sidebar: product1_sidebar
permalink: sw_cff.html
folder: swfiles
mathjax: true
---
  returns the atomic charge form factor values for X-ray scattering
 
  [formFactVal, coeff] = SW_CFF(atomName, {Q})
 
  The provided form factor values at Q=0 are normalized to Z.
 
  Input:
 
  atomName      String, contains the name of the ion in using the symbol of
                the element following the charge, e.g. Cr3+. It can be also
                the coefficients to calculate f. If the string contains
                whitespace, the first word will be used as input.
  Q             Momentum transfer in Angstrom^-1 units with dimensions of
                [1 nQ] or [3 nQ], optional.
 
  Output:
 
  formFactVal   Value of form factor, evaluated at the Q points if Q is
                defined.
  coeff         Form factor coefficients according to the following
                formula:
                    f0(Qs) = c + SUM a_i*EXP(-b_i*(Qs^2))
                                 i=1,5
                where Qs = Q/(4*pi) and [a_1 b_1 a_2 b_2 ... c] are the
                coefficients.
 
  See also sw_mff.
 
