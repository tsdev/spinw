---
title: base( )
keywords: sample
summary: "sets the basis vectors of an swplot figure"
sidebar: product1_sidebar
permalink: base.html
folder: +swplot
mathjax: true
---
  sets the basis vectors of an swplot figure
 
  SWPLOT.BASE(BV, {hFigure})
 
  BV is a matrix with dimensions of [3 3] and contains the three basis
  vectors of the new coordinate system as column vectors.
 
  SWPLOT.BASE(obj, {hFigure})
 
  obj is a spinw object that defines the swplot coordinate system as
  lattice units.
 
  BV = SWPLOT.BASE
 
  Returns the basis vectors stored in the swplot figure.
 
  Input:
 
  BV            Either a 3x3 matrix of the new basis vectors or a spinw
                object where the new basis vectors will be the lattice
                units of the stored crystal.
  hFigure       Handle of the swplot figure. Default is the selected
                figure.
 
  See also SWPLOT.PLOT.
 
