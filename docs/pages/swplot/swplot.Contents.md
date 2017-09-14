---
{title: swplot package, link: swplot package, summary: package for 3D plotting, keywords: sample,
  sidebar: sw_sidebar, permalink: swplot.html, folder: swplot, mathjax: 'true'}

---
 
The library contains functions that can create and control plot of
3D objects. It contains low level function to plot specific shapes on any
figure (cylinder, line, arrow, ellipsoid, text) that are vectorized.
Moreover it provides the plot command to plot multiple objects with
specified color, transparency, label and legend on the swplot figure that
has rotation and translations built in that are better than the default
3D rotation with mouse in Matlab. Moreover it has high level commands to
plot crystals from SpinW objects (plotatom, plotmag, plotion, plotbond,
plotbase, plotcell, plotchem).
 
### Files
 
#### Basic 3D objects
 
* [swplot.arrow](swplot_arrow.html) draws a 3D arrow using patch
* [swplot.circle](swplot_circle.html) creates a circle surface in 3 dimensions
* [swplot.cylinder](swplot_cylinder.html) draws a closed/open 3D cylinder
* [swplot.ellipsoid](swplot_ellipsoid.html) draw ellipsoid
* [swplot.line](swplot_line.html) draws a 3D line using patch
* [swplot.polyhedron](swplot_polyhedron.html) draw convex polyhedra or polynom from vertex list
* [swplot.text](swplot_text.html) draws a text at a point in 3D
 
#### Creating and modifying swplot figure
 
* [swplot.activefigure](swplot_activefigure.html) returns the handle of the active swplot figure
* [swplot.add](swplot_add.html) adds a graphical object to the hgtransform of an swplot figure
* [swplot.base](swplot_base.html) sets the basis vectors of an swplot figure
* [swplot.clear](swplot_clear.html) clear swplot figure
* [swplot.close](swplot_close.html) close swplot figure
* [swplot.delete](swplot_delete.html) deletes objects and their data on swplot figure
* [swplot.export](swplot_export.html) exports swplot figure into raster/vector image
* [swplot.figure](swplot_figure.html) creates swplot figure
* [swplot.findobj](swplot_findobj.html) finds object data on swplot figure
* [swplot.getdata](swplot_getdata.html) gets the data stored in an swplot figure
* [swplot.ishg](swplot_ishg.html) does the swplot figure uses hgtransform
* [swplot.legend](swplot_legend.html) draws legend to the swplot figure
* [swplot.mouse](swplot_mouse.html) adds mouse callbacks to swplot figure
* [swplot.plot](swplot_plot.html) plots objects to swplot figure
* [swplot.tooltip](swplot_tooltip.html) creates tooltip axis on swplot figure
* [swplot.transform](swplot_transform.html) transform objects on swplot figure
* [swplot.translate](swplot_translate.html) translate objects on swplot figure
* [swplot.view](swplot_view.html) control 3D view of swplot
* [swplot.zoom](swplot_zoom.html) zooming objects on swplot figure
 
#### Plotting SpinW object on swplot figure
 
* [swplot.plotatom](swplot_plotatom.html) plots crystal structure
* [swplot.plotbase](swplot_plotbase.html) plots the edges of unit cells on swplot figure
* [swplot.plotbond](swplot_plotbond.html) plots magnetic bonds
* [swplot.plotcell](swplot_plotcell.html) plots the edges of unit cells on swplot figure
* [swplot.plotchem](swplot_plotchem.html) plots polyhedra or chemical bonds
* [swplot.plotion](swplot_plotion.html) plots magnetic ion properties
* [swplot.plotmag](swplot_plotmag.html) plots magnetic structure
 
#### Related functions
 
* [swplot.color](swplot_color.html) generates RGB code from color name string
* [swplot.icomesh](swplot_icomesh.html) creates mesh by subdividing icosahedron faces
* [swplot.logo](swplot_logo.html) creates the logo of SpinW and saves to a .png file
* [swplot.patchfacefcn](swplot_patchfacefcn.html) callback function for patch face selection
* [swplot.raytriangle](swplot_raytriangle.html) finds if a ray crosses a triangle
* [swplot.setrangegui](swplot_setrangegui.html) shows a window to change the plotting range of an swplot figure
* [swplot.subfigure](swplot_subfigure.html) changes position of figure window on the screen
* [swplot.tooltipcallback](swplot_tooltipcallback.html) call tooltip on clicking an object on an swplot figure
* [swplot.tooltipstring](swplot_tooltipstring.html) generate tooltip string from the data of a graphical object
 

