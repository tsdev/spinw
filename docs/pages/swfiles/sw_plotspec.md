---
{title: sw_plotspec( ), summary: plots spin wave spectrum, keywords: sample, sidebar: sw_sidebar,
  permalink: sw_plotspec.html, folder: swfiles, mathjax: 'true'}

---
plots spin wave spectrum
 
[fHandle, pHandle] = SW_PLOTSPEC(spectra, 'option1', value1 ...)
 
The color plots using cm_inferno() as a default colormap, change the
default colormap using [swpref.setpref](swpref_setpref.html)('colormap',@my_colomap) command.
 
Options:
 
mode      Choose the type of plot, either a string (or number):
              'disp'  dispersion (1),
              'int'   intensity of the correlation functions (2),
              'color' convoluted spectrum (3),
              'fancy' FANCY PLOT MODE (default, 4).
imag      Whether to plot the imaginary values of the dispersion
          and the correlation functions. For convoluted spectra, if true,
          the imaginary part is plotted. Default is false.
aHandle   Handle of the axis object for plotting, if undefined the
          active axis will be used (gca).
colorbar  Plot colorbar for dispersion and intensity, default is true.
nCol      Number of colors in the colormap, default is 500.
dashed    Whether to plot dashed vertical line between multiple linear
          scans. Defult is false.
dE        FWHM value of convoluted Gaussian in energy to simulate finite
          energy resolution. Only works for mode=3. If zero, no
          convolution performed. Default is 0.
fontSize  Font size on the plot, default is 14 pt.
colormap  Colormap for plotting, default is stored in 
          [swpref.getpref](swpref_getpref.html)('colormap'). For single plot and for multiple
          plot it will be a continuous scale from white to different
          color. This is the 'auto' mode. Also colormap can be given
          directly using standard colormaps, like @jet. To overplot
          multiple spectras 'colormap' option will be a matrix, with
          dimensions [3 nConv], where every column defines a color for
          the maximum intensity. It is also used for plotting dispersion
          curves. In case a single color all dispersion curves have the
          same color (e.g. [255 0 0] for red), or as many colors as
          dispersion curves (dimensions are [3 nMode]), or any colormap
          can be given, like @jet. In this case every mode will have
          different colors, the color is determined from the index of the
          mode. Default is 'auto'.
sortMode  Sorting the modes before plotting. Default is false.
axLim     Upper limit for y axis (mode 1,2) or z axis (mode 3), default
          is 'auto'. For color plot of multiple cross section the c axis
          cannot be changed after the plot.
legend    Whether to plot legend for multiple convoluted spectras,
          default is true.
title     Whether to plot figure title, default is true.
twin      Select which twins to plot for omega plots, default plots all
          twins, dimensions are [1 nTwinToPlot].
lineStyle Line style for line plots (dispersion and intensity), default
          is {'-' 'o-' '--'}. For example '--' gives dashed lines.
lineWidth Line width of line plots, default is 0.5 point.
log       Plot 10based logarithmic intensity, default is false.
plotf     Plot function for color plot. Default is @surf.
maxPatch  Maximum number of pixels that can be plotted using the patch()
          function within sw_surf(). Using patch for color plot can be
          slow on older machines, but the figure can be exported
          afterwards as a vector graphics, using the print() function.
          Default is 1000.
norm      If true, the convolution with a Gaussian function (in case of
          non-zero 'dE' option) keeps the energy integrated intensity. If
          false the amplitude is kept constant. Default is the input
          spectra.norm value.
x0        Vector with two numbers [x0_min x0_max]. By default the x range
          of the plots is [0 1] irrespective of the given Q points. To
          change this the lower and upper limits can be given here.
qlabel    Provide a list of strings for the Q points between linear
          segments.
dat       Experimental data points to plot over the calculated spectrum.
          Can be either the name of a data file that contain the
          experimentally fitted dispersion (needs to have the same format
          as the input for [spinw](spinw.html).fitspec()), or it is a structure that
          contains the already imported data using [sw_readtable()](sw_readtable.html), for
          example:
              T = [sw_readtable](sw_readtable.html)('myExpData.txt','	');
              [sw_plotspec](sw_plotspec.html)(spectra,'dat',T);
ddat      Maximum distance between any Q point in the simulated spectrum
          and an experimental data point in A^-1. If an experimental data
          point is further from any Q point of the simulated spectrum
          thant ddat, it won't be plotted. Default value is 0.01 A^-1.
 
Output:
 
fHandle   Handle of the plot figure.
pHandle   Handle of the graphics objects on the figure.
 
See also SPINW.PLOT, SPINW.SPINWAVE, SW_SURF, SW_LABEL.
 

