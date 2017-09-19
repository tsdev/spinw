---
{title: sw_multicolor( ), link: sw_multicolor, summary: creates RGB color data for
    multiple 2D overlapping plots, keywords: sample, sidebar: sw_sidebar, permalink: sw_multicolor.html,
  folder: swfiles, mathjax: 'true'}

---

### Syntax

` `

### Description



### Examples

Plotting of two random matrices (dimensions are [100 100]) with
red and blue colors:
    cMat = sw_multicolor(rand(100,100,2),[1 0;0 1;0 0],[0 1]);
    image(cMat);

### Input Arguments

% `vMat`
:  Matrix that contains the input 2D data, dimensions are

% `[d1`
:Plot]. Where each plot has a dimensions of [d1 d2].

% `cMap`
:  Cell of colormap functions containing used for the different

% `overlayed`
:d plots. For example:

% `cMap`
:@copper @gray}.

% `cLim`
:  Maximum and minimum values of the color map. Values in vMat

% `smaller`
:than the minimum and larger than the maximum will be

% `shown`
:th the minimum and maximum values in the colormap

% `respectively.`
:vely. The dimensions of cLim is [1 2].

% `nCol`
:  Number of colors in the colormap. Optional, default value is

% ``
:

% `flipud`
:  If true the colormaps are inverted. Optional, default is false.

### Output Arguments

cMat      Matrix with equal dimensions to the input times three for the
red, green and blue channels, dimensions are [d1 d2 3].

