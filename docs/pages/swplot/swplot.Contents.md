---
{title: swplot package, link: swplot package, summary: package for 3D plotting, keywords: sample,
  sidebar: sw_sidebar, permalink: swplot, folder: swplot, mathjax: 'true'}

---
 
The package contains functions that can create and control plot of 3D
objects. It contains low level function to plot specific shapes on any
figure (cylinder, line, arrow, ellipsoid, text) that are vectorized.
Moreover it provides the plot command to plot multiple objects with
specified color, transparency, label and legend on the swplot figure that
has rotation and translations built in that are better than the default
3D rotation with mouse in Matlab. Moreover it has high level commands to
plot crystals from SpinW objects (plotatom, plotmag, plotion, plotbond,
plotbase, plotcell, plotchem).
 
### Files
 
#### Basic 3D objects
 
* [swplot.arrow](swplot_arrow) creates a 3D arrow patch
* [swplot.circle](swplot_circle) creates a 3D circle surface patch
* [swplot.cylinder](swplot_cylinder) creates a closed/open 3D cylinder patch
* [swplot.ellipsoid](swplot_ellipsoid) creates a 3D ellipsoid patch
* [swplot.line](swplot_line) creates 3D line patch
* [swplot.polyhedron](swplot_polyhedron) draw convex polyhedra or polynom from vertex list
* [swplot.text](swplot_text) draws a text at a point in 3D
 
#### Creating and modifying swplot figure
 
* [swplot.activefigure](swplot_activefigure) returns the handle of the active swplot figure
* [swplot.add](swplot_add) adds a graphical object to an swplot figure
* [swplot.base](swplot_base) sets the basis vectors of an swplot figure
* [swplot.clear](swplot_clear) clears swplot figure
* [swplot.close](swplot_close) closes swplot figure
* [swplot.delete](swplot_delete) deletes objects and corresponding data from swplot figure
* [swplot.export](swplot_export) exports swplot figure into raster/vector image
* [swplot.figure](swplot_figure) creates swplot figure
* [swplot.findobj](swplot_findobj) finds object data on swplot figure
* [swplot.getdata](swplot_getdata) gets the data stored in an swplot figure
* [swplot.ishg](swplot_ishg) does the swplot figure uses hgtransform
* [swplot.legend](swplot_legend) adds legend to the swplot figure
* [swplot.mouse](swplot_mouse) adds mouse callbacks to swplot figure
* [swplot.plot](swplot_plot) plots objects to swplot figure
* [swplot.tooltip](swplot_tooltip) creates tooltip axis on swplot figure
* [swplot.transform](swplot_transform) transform objects on swplot figure
* [swplot.translate](swplot_translate) translate objects on swplot figure
* [swplot.view](swplot_view) control 3D view of swplot
* [swplot.zoom](swplot_zoom) zooming objects on swplot figure
 
#### Plotting SpinW object on swplot figure
 
* [swplot.plotatom](swplot_plotatom) plots crystal structure
* [swplot.plotbase](swplot_plotbase) plots the edges of unit cells on swplot figure
* [swplot.plotbond](swplot_plotbond) plots magnetic bonds
* [swplot.plotcell](swplot_plotcell) plots the edges of unit cells on swplot figure
* [swplot.plotchem](swplot_plotchem) plots polyhedra or chemical bonds
* [swplot.plotion](swplot_plotion) plots magnetic ion properties
* [swplot.plotmag](swplot_plotmag) plots magnetic structure
 
#### Related functions
 
* [swplot.color](swplot_color) generates RGB code from color name
* [swplot.icomesh](swplot_icomesh) creates mesh by subdividing icosahedron faces
* [swplot.logo](swplot_logo) creates the SpinW logo
* [swplot.patchfacefcn](swplot_patchfacefcn) callback function for patch face selection
* [swplot.raytriangle](swplot_raytriangle) finds if a ray crosses a triangle
* [swplot.setrangegui](swplot_setrangegui) shows a window to change the plotting range of an swplot figure
* [swplot.subfigure](swplot_subfigure) changes position of figure window on the screen
* [swplot.subplot](swplot_subplot) create subplots with variable gaps between axes
* [swplot.tooltipcallback](swplot_tooltipcallback) call tooltip on clicking an object on an swplot figure
* [swplot.tooltipstring](swplot_tooltipstring) generate tooltip string from the data of a graphical object
 

{% include links.html %}
