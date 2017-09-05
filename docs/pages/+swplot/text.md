---
title: text( )
keywords: sample
summary: "draws a text at a point in 3D"
sidebar: product1_sidebar
permalink: text.html
folder: +swplot
mathjax: true
---
  draws a text at a point in 3D
 
  hText = SWPLOT.TEXT(r, string, {fontSize})
 
  hText = SWPLOT.TEXT(handle,...)
 
  Handle of an axes object that selects an axis to plot.
 
  Input:
 
  handle    Handle of an axis object.
  r         Coordinate of the center of the text for a single text or
            matrix with dimensions [3 nText] for multiple text.
  string    String contains the text or cell of strings to plot multiple
            text.
  fontSize  Font size in pt, default value is stored in
            swpref.getpref('fontsize')
 
  See also TEXT.
 
