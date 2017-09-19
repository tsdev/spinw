---
{title: swplot.plotcell, link: swplot.plotcell, summary: plots the edges of unit cells
    on swplot figure, keywords: sample, sidebar: sw_sidebar, permalink: swplot_plotcell.html,
  folder: swplot, mathjax: 'true'}

---

### Syntax

` `

### Description

 

### Name-Value Pair Arguments

% `range`
:  Plotting range of the lattice parameters in lattice units,

% `dimensions`
:ns are [3 2]. For example to plot the first unit cell,

% `use:`
:1;0 1;0 1]. Also the number unit cells can be given

% `along`
:e a, b and c directions: [2 1 2], that is equivalent to

% `[0`
:;0 2]. Default is the single unit cell.

% `mode`
:  String, determined how the cells are plotted:
 gle'    A single unit cell is plotted at the origin.
 ide'    Unit cells are plotted inside the given
         range. Default.
 side'   Unit cells are plotted inclusive the given
             range.

% `figure`
:  Handle of the swplot figure. Default is the selected figure.

% `color`
:  Color of the lines as a string or row vector with 3 elements, 

% `default`
:value is black ('auto').

% `lineStyle`
:e Line style of the cell, default is '--'.

% `lineWdith`
:h Line width of the cell, default is 1.

% `translate`
:e If true, all plot objects will be translated to the figure

% `center.`
:Default is false.

% `zoom`
:  If true, figure will be automatically zoomed to the ideal size.

% `Default`
:is false.

