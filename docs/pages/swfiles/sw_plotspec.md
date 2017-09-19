---
{title: sw_plotspec( ), link: sw_plotspec, summary: plots spin wave spectrum, keywords: sample,
  sidebar: sw_sidebar, permalink: sw_plotspec.html, folder: swfiles, mathjax: 'true'}

---

### Syntax

` `

### Description

default colormap using swpref.setpref('colormap',@my_colomap) command.
 

### Name-Value Pair Arguments

% `mode`
:  Choose the type of plot, either a string (or number):
 p'  dispersion (1),
 '   intensity of the correlation functions (2),
 or' convoluted spectrum (3),
 cy' FANCY PLOT MODE (default, 4).

% `imag`
:  Whether to plot the imaginary values of the dispersion

% `and`
:correlation functions. For convoluted spectra, if true,

% `the`
:inary part is plotted. Default is false.

% `aHandle`
:  Handle of the axis object for plotting, if undefined the

% `active`
:xis will be used (gca).

% `colorbar`
:  Plot colorbar for dispersion and intensity, default is true.

% `nCol`
:  Number of colors in the colormap, default is 500.

% `dashed`
:  Whether to plot dashed vertical line between multiple linear

% `scans.`
:efult is false.

% `dE`
:  FWHM value of convoluted Gaussian in energy to simulate finite

% `energy`
:esolution. Only works for mode=3. If zero, no

% `convolution`
:ion performed. Default is 0.

% `fontSize`
:  Font size on the plot, default is 14 pt.

% `colormap`
:  Colormap for plotting, default is stored in 

% `swpref.getpref('colormap').`
:etpref('colormap'). For single plot and for multiple

% `plot`
:will be a continuous scale from white to different

% `color.`
:his is the 'auto' mode. Also colormap can be given

% `directly`
: using standard colormaps, like @jet. To overplot

% `multiple`
: spectras 'colormap' option will be a matrix, with

% `dimensions`
:ns [3 nConv], where every column defines a color for

% `the`
:mum intensity. It is also used for plotting dispersion

% `curves.`
:In case a single color all dispersion curves have the

% `same`
:or (e.g. [255 0 0] for red), or as many colors as

% `dispersion`
:on curves (dimensions are [3 nMode]), or any colormap

% `can`
:iven, like @jet. In this case every mode will have

% `different`
:t colors, the color is determined from the index of the

% `mode.`
:fault is 'auto'.

% `sortMode`
:  Sorting the modes before plotting. Default is false.

% `axLim`
:  Upper limit for y axis (mode 1,2) or z axis (mode 3), default

% `is`
:'. For color plot of multiple cross section the c axis

% `cannot`
:e changed after the plot.

% `legend`
:  Whether to plot legend for multiple convoluted spectras,

% `default`
:is true.

% `title`
:  Whether to plot figure title, default is true.

% `twin`
:  Select which twins to plot for omega plots, default plots all

% `twins,`
:imensions are [1 nTwinToPlot].

% `lineStyle`
:e Line style for line plots (dispersion and intensity), default

% `is`
:'o-' '--'}. For example '--' gives dashed lines.

% `lineWidth`
:h Line width of line plots, default is 0.5 point.

% `log`
:  Plot 10based logarithmic intensity, default is false.

% `plotf`
:  Plot function for color plot. Default is @surf.

% `maxPatch`
:  Maximum number of pixels that can be plotted using the patch()

% `function`
: within sw_surf(). Using patch for color plot can be

% `slow`
:older machines, but the figure can be exported

% `afterwards`
:ds as a vector graphics, using the print() function.

% `Default`
:is 1000.

% `norm`
:  If true, the convolution with a Gaussian function (in case of

% `non-zero`
: 'dE' option) keeps the energy integrated intensity. If

% `false`
:e amplitude is kept constant. Default is the input

% `spectra.norm`
:norm value.

% `x0`
:  Vector with two numbers [x0_min x0_max]. By default the x range

% `of`
:lots is [0 1] irrespective of the given Q points. To

% `change`
:his the lower and upper limits can be given here.

% `qlabel`
:  Provide a list of strings for the Q points between linear

% ``
:.

% `dat`
:  Experimental data points to plot over the calculated spectrum.

% `Can`
:ither the name of a data file that contain the

% `experimentally`
:ntally fitted dispersion (needs to have the same format

% `as`
:nput for spinw.fitspec()), or it is a structure that

% `contains`
: the already imported data using sw_readtable(), for

% ``
:
 sw_readtable('myExpData.txt','\t');
 lotspec(spectra,'dat',T);

% `ddat`
:  Maximum distance between any Q point in the simulated spectrum

% `and`
:xperimental data point in A^-1. If an experimental data

% `point`
: further from any Q point of the simulated spectrum

% `thant`
:at, it won't be plotted. Default value is 0.01 A^-1.

### Output Arguments

fHandle   Handle of the plot figure.
pHandle   Handle of the graphics objects on the figure.

### See Also

[spinw.plot](spinw_plot.html), [spinw.spinwave](spinw_spinwave.html), [sw_surf] and [sw_label](sw_label.html)

