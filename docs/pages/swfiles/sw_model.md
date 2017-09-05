---
title: sw_model( )
keywords: sample
summary: "creates different predefined spin models"
sidebar: product1_sidebar
permalink: sw_model.html
folder: swfiles
mathjax: true
---
  creates different predefined spin models
 
  obj = SW_MODEL(model, param, {fid})
 
  Input:
 
  model     String, name of the model, one of the following:
                'triAF'     Triangular lattice Heisenberg antiferromagnet
                            in the ab plane (a=b=3 Angstrom), with gamma =
                            120 deg angle and optimised magnetic structure.
                            Arbitrary number of Heisenberg interaction can
                            be defined, param(1) gives the value of 1st
                            neighbor interaction, param(2) the second etc.
                'squareAF'  Square lattice antiferromagnet.
                'chain'     Chain with further neighbor interactions.
 
  param     Input parameters of the model, depending on which is selected.
  fid       Where to print the text output. Default is 1 to print to the
            Command Window. Optional.
 
  Output:
 
  obj       spinw class object with the selected model.
 
  See also SPINW.
 
