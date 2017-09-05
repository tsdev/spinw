---
title: tooltip( )
keywords: sample
summary: "creates tooltip axis on swplot figure"
sidebar: product1_sidebar
permalink: tooltip.html
folder: +swplot
mathjax: true
---
  creates tooltip axis on swplot figure
 
  SWPLOT.TOOLTIP({text}, {hFigure}, {window})
 
  status = SWPLOT.TOOLTIP
 
  Input:
 
  text          String, if it is 'on'/'off' the tooltip will be switched
                on/off. Otherwise the text will be shown in the tooltip.
                Default is 'on'.
  hFigure       Handle of the swplot figure. Default is the selected
                figure.
  window        If true, the tooltips will be shown in a separate window.
                Default is false.
 
  Output:
 
  status        String, can be 'on' or 'off'.
 
  Example:
    swplot.figure
    swplot.addcircle([0 0 0],[0 0 1],1)
    swplot.tooltip
 
