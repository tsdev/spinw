---
{title: sw_plotspec, link: sw_plotspec, summary: plots spectrum, keywords: sample,
  sidebar: sw_sidebar, permalink: sw_plotspec, folder: swfiles, mathjax: 'true'}

---
  
### Syntax
  
`[fhandle, phandle] = sw_plotspec(spectra,Name,Value)`
  
### Description
  
`[fhandle, phandle] = sw_plotspec(spectra,Name,Value)` plots excitation
spectrum that is calculated either by [spinw.spinwave](spinw_spinwave) or
[spinw.powspec](spinw_powspec). It can plot dispersion or intensities as line plots or
the energy binned spectrum as a color plot. The color plots uses
[cm_inferno] as a default colormap. To change the default colormap use
the `swpref.setpref('colormap',@my_colomap)` command. The function is
able to plot the spectrum if it is calculated along a path in the
Brillouin-zone and display the labels of the high symmetry Brillouon-zone
points.
  
### Name-Value Pair Arguments
  
`'mode'`
: Choose the type of plot using the following strings:
  * `'disp'`  Plot dispersion as line plot.
  * `'int'`   PLot intensity of each mode as line plot.
  * `'color'` Color plot of energy binned spectrum.
  * `'auto'`  Auto plot mode that tries to determine the best
              parameteres, default.
  
`'imag'`
: If `true` also the imaginary part of the dispersion
  and the correlation function values will be shown as red lines on top
  the real values. For color plot if `true` only the imaginary part of
  the binned data will be shown. Default value is `false`.
  
`'aHandle'`
: Handle of the axis object which will show the plot. If undefined the
  active axis will be used, see [gca].
  
`'colorbar'`
: Plot colorbar for dispersion and intensity, default value is `true`.
  
`'nCol'`
: Number of colors in the colormap, default value is 500.
  
`'dashed'`
: If `true` dashed vertical lines between linear $$Q$$ segments will be
  shown. Default is `false`.
  
`'dE'`
: If given, a Gaussian will be convoluted with the binned data to simulate finite
  energy resolution. Only works if `mode=3`. If zero, no convolution
  performed. Default value is 0.
  
`'fontSize'`
: Font size in pt for the labels on the plot, default value is 14 pt.
  
`'colormap'`
: Colormap for plotting, default value is stored in 
  `swpref.getpref('colormap')`. For single plot and for multiple plot it
  will be a continuous scale from white to different color. This is the
  `'auto'` mode. Also colormap can be given directly using standard
  colormaps as function handles, e.g. `@jet`. To overplot multiple
  spectra `colormap` option will be a matrix, with dimensions [3 nConv],
  where every column defines a color for the maximum intensity. It is
  also used for plotting dispersion curves. In case a single color all
  dispersion curves have the same color (e.g. `[255 0 0]` for red), or as
  many colors as dispersion curves in a matrix with dimensions of
  $$[3\times n_{mode}]$$ or as a colormap function handle. In this case
  every mode will have different color and the color is determined from
  the index of the mode after the colormap is applied. Default value is
  `'auto'`.
  
`'sortMode'`
: Sorting the modes before plotting. Default is `false`. Can improve the
  quality of the dispersion line plots if modes are crossing.
  
`'axLim'`
: Upper limit for energy axis (for `mode` 1,2) or color axis (for `mode`
  3), default value is `'auto'`. For color plot of multiple spectra
  the color axis cannot be changed after the plot.
  
`'legend'`
: Whether to plot legend for multiple convoluted spectras,
  default value is `true`.
  
`'title'`
: If `true` a title will be added to the figure, default value is `true`.
  
`'twin'`
: Select which twins to be plotted for dispersion plots, by default the
  spectrum corresponding to all twins will be plotted. The dimensions are
  $$[1\times n_{twinToPlot}]$$.
  
`'lineStyle'`
: Line style for line plots (dispersion and intensity), default value
  `{'-' 'o-' '--'}` for plotting modes that correspond to line style of
  $$S(Q,\omega)$$, $$S(Q+k,\omega)$$ and $$S(Q-k,\omega)$$ cross modes in case
  of incommensurate magnetic systems. For commensurate systems only thte
  first string in the cell will be considered. For example '--' gives
  dashed lines.
  
`'lineWidth'`
: Line width of line plots, default value is 0.5 pt.
  
`'log'`
: If true, the 10-based logarithmic intensity will be plotted, default
  value is `false`.
  
`'plotf'`
: Function handle of the plot function for color plot. Default is
  `surf`.
  
`'maxPatch'`
: Maximum number of pixels that can be plotted using the [patch]
  function within [sw_surf]. Using [patch] for color plot can be
  slow on older machines, but the figure can be exported
  afterwards as a vector graphics, using the [print] function.
  Default value is 1000.
  
`'norm'`
: If true, the convolution with a Gaussian function (in case of
  non-zero `dE` parameter) keeps the energy integrated intensity. If
  `false` the amplitude is kept constant. Default is determined by the
  value stored in the input `spectra.norm`.
  
`'x0'`
: Row vector with two numbers `[x0_min x0_max]`. By default the $$x$$ range
  of the plot is `[0 1]` irrespective of the values of the $$Q$$ values. To
  change this the lower and upper limits can be given here.
  
`'qlabel'`
: Provide a list of strings for the special $$Q$$ points along the path in
  the Brillouin zone, e.g. `{'X' '\Gamma' 'M' 'K' '\Gamma'}`.
  
`'dat'`
: Experimental data points to plot over the calculated spectrum.
  Can be either the name of a data file that contains the
  experimentally fitted dispersion (needs to have the same format
  as the input file of [spinw.fitspec](spinw_fitspec) see help for details on the file
  format), or it is a structure that contains the already imported data
  using the [sw_readtable](sw_readtable) function, e.g.
 
  ```matlab
  T = sw_readtable('myExpData.txt','\t');
  sw_plotspec(spectra,'dat',T);
  ```
  
`'ddat'`
: Maximum distance between any $$Q$$ point in the simulated spectrum
  and an experimental data point in Å$$^{-1}$$ unit. If an
  experimental data point is further from any $$Q$$ point than the given 
  limit, it will be omitted. Default value is 0.01.
  
### Output Arguments
  
`fHandle`
: Handle of the plot figure.
 
`pHandle`
: Vector that contains the handle of the graphics objects on the figure.
  
### See Also
  
[spinw.plot](spinw_plot) \| [spinw.spinwave](spinw_spinwave) \| [sw_surf] \| [sw_label]
 
*[FWHM]: Full Width at Half Maximum
 

{% include links.html %}
