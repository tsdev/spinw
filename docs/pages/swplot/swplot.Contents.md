---
{title: Package swplot, summary: The SWPLOT library contains functions that can create
    and control plot of, keywords: sample, sidebar: sw_sidebar, permalink: swplot.html,
  folder: swplot, mathjax: 'true'}

---
The SWPLOT library contains functions that can create and control plot of
3D objects. It contains low level function to plot specific shapes on any
figure (cylinder, line, arrow, ellipsoid, text) that are vectorized.
Moreover it provides the plot command to plot multiple objects with
specified color, transparency, label and legend on the [swplot figure](swplot_figure.html) that
has rotation and translations built in that are better than the default
3D rotation with mouse in Matlab. Moreover it has high level commands to
plot crystals from SpinW objects (plotatom, plotmag, plotion, plotbond,
plotbase, plotcell, plotchem).
 
### Files
 
#### Basic 3D objects
 

* [[swplot.arrow](swplot_arrow.html)](/swplot_arrow) draws a 3D arrow using patch
* [[swplot.circle](swplot_circle.html)](/swplot_circle) creates a circle surface in 3 dimensions
* [[swplot.cylinder](swplot_cylinder.html)](/swplot_cylinder) draws a closed/open 3D cylinder
* [[swplot.ellipsoid](swplot_ellipsoid.html)](/swplot_ellipsoid) draw ellipsoid
* [[swplot.line](swplot_line.html)](/swplot_line) draws a 3D line using patch
* [[swplot.polyhedron](swplot_polyhedron.html)](/swplot_polyhedron) draw convex polyhedra or polynom from vertex list
* [[swplot.text](swplot_text.html)](/swplot_text) draws a text at a point in 3D
 
#### Creating and modifying [swplot figure](swplot_figure.html)
 

* [[swplot.activefigure](swplot_activefigure.html)](/swplot_activefigure) returns the handle of the active [swplot figure](swplot_figure.html)
* [[swplot.add](swplot_add.html)](/swplot_add) adds a graphical object to the hgtransform of an [swplot figure](swplot_figure.html)
* [[swplot.base](swplot_base.html)](/swplot_base) sets the basis vectors of an [swplot figure](swplot_figure.html)
* [[swplot.clear](swplot_clear.html)](/swplot_clear) clear [swplot figure](swplot_figure.html)
* [[swplot.close](swplot_close.html)](/swplot_close) close [swplot figure](swplot_figure.html)
* [[swplot.delete](swplot_delete.html)](/swplot_delete) deletes objects and their data on [swplot figure](swplot_figure.html)
* [[swplot.export](swplot_export.html)](/swplot_export) exports [swplot figure](swplot_figure.html) into raster/vector image
* [[swplot.figure](swplot_figure.html)](/swplot_figure) creates [swplot figure](swplot_figure.html)
* [[swplot.findobj](swplot_findobj.html)](/swplot_findobj) finds object data on [swplot figure](swplot_figure.html)
* [[swplot.getdata](swplot_getdata.html)](/swplot_getdata) gets the data stored in an [swplot figure](swplot_figure.html)
* [[swplot.ishg](swplot_ishg.html)](/swplot_ishg) does the [swplot figure](swplot_figure.html) uses hgtransform
* [[swplot.legend](swplot_legend.html)](/swplot_legend) draws legend to the [swplot figure](swplot_figure.html)
* [[swplot.mouse](swplot_mouse.html)](/swplot_mouse) adds mouse callbacks to [swplot figure](swplot_figure.html)
* [[swplot.plot](swplot_plot.html)](/swplot_plot) plots objects to [swplot figure](swplot_figure.html)
* [[swplot.tooltip](swplot_tooltip.html)](/swplot_tooltip) creates tooltip axis on [swplot figure](swplot_figure.html)
* [[swplot.transform](swplot_transform.html)](/swplot_transform) transform objects on [swplot figure](swplot_figure.html)
* [[swplot.translate](swplot_translate.html)](/swplot_translate) translate objects on [swplot figure](swplot_figure.html)
* [[swplot.view](swplot_view.html)](/swplot_view) control 3D view of swplot
* [[swplot.zoom](swplot_zoom.html)](/swplot_zoom) zooming objects on [swplot figure](swplot_figure.html)
 
#### Plotting SpinW object on [swplot figure](swplot_figure.html)
 

* [[swplot.plotatom](swplot_plotatom.html)](/swplot_plotatom) plots crystal structure
* [[swplot.plotbase](swplot_plotbase.html)](/swplot_plotbase) plots the edges of unit cells on [swplot figure](swplot_figure.html)
* [[swplot.plotbond](swplot_plotbond.html)](/swplot_plotbond) plots magnetic bonds
* [[swplot.plotcell](swplot_plotcell.html)](/swplot_plotcell) plots the edges of unit cells on [swplot figure](swplot_figure.html)
* [[swplot.plotchem](swplot_plotchem.html)](/swplot_plotchem) plots polyhedra or chemical bonds
* [[swplot.plotion](swplot_plotion.html)](/swplot_plotion) plots magnetic ion properties
* [[swplot.plotmag](swplot_plotmag.html)](/swplot_plotmag) plots magnetic structure
 
#### Related functions
 

* [[swplot.color](swplot_color.html)](/swplot_color) generates RGB code from color name string
* [[swplot.icomesh](swplot_icomesh.html)](/swplot_icomesh) creates mesh by subdividing icosahedron faces
* [[swplot.logo](swplot_logo.html)](/swplot_logo) creates the logo of SpinW and saves to a .png file
* [[swplot.patchfacefcn](swplot_patchfacefcn.html)](/swplot_patchfacefcn) callback function for patch face selection
* [[swplot.raytriangle](swplot_raytriangle.html)](/swplot_raytriangle) finds if a ray crosses a triangle
* [[swplot.setrangegui](swplot_setrangegui.html)](/swplot_setrangegui) shows a window to change the plotting range of an [swplot figure](swplot_figure.html)
* [[swplot.subfigure](swplot_subfigure.html)](/swplot_subfigure) changes position of figure window on the screen
* [[swplot.tooltipcallback](swplot_tooltipcallback.html)](/swplot_tooltipcallback) call tooltip on clicking an object on an [swplot figure](swplot_figure.html)
* [[swplot.tooltipstring](swplot_tooltipstring.html)](/swplot_tooltipstring) generate tooltip string from the data of a graphical object
 

