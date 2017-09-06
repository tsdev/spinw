---
{title: Contents( ), keywords: sample, summary: The SWPLOT library contains functions that can create and control plot of,
  sidebar: sw_sidebar, permalink: +swplot_Contents.html, folder: +swplot, mathjax: 'true'}

---
  The SWPLOT library contains functions that can create and control plot of
  3D objects. It contains low level function to plot specific shapes on any
  figure (cylinder, line, arrow, ellipsoid, text) that are vectorized.
  Moreover it provides the plot command to plot multiple objects with
  specified color, transparency, label and legend on the swplot figure that
  has rotation and translations built in that are better than the default
  3D rotation with mouse in Matlab. Moreover it has high level commands to
  plot crystals from SpinW objects (plotatom, plotmag, plotion, plotbond,
  plotbase, plotcell, plotchem).
 
  Files
 
  Vectorized plotting of basic 3D geometrical objects using the patch
  command and triangular patches:
 
    arrow           - draws a 3D arrow using patch
    circle          - creates a circle surface in 3 dimensions
    cylinder        - draws a closed/open 3D cylinder
    ellipsoid       - draw ellipsoid
    line            - draws a 3D line using patch
    polyhedron      - draw convex polyhedra or polynom from vertex list
    text            - draws a text at a point in 3D
 
  Creating and modifying swplot figure:
 
    activefigure    - returns the handle of the active swplot figure
    add             - adds a graphical object to the hgtransform of an swplot figure
    base            - sets the basis vectors of an swplot figure
    clear           - clear swplot figure
    close           - close swplot figure
    delete          - deletes objects and their data on swplot figure
    export          - exports swplot figure into raster/vector image
    figure          - creates swplot figure
    findobj         - finds object data on swplot figure
    getdata         - gets the data stored in an swplot figure
    ishg            - does the swplot figure uses hgtransform
    legend          - draws legend to the swplot figure
    mouse           - adds mouse callbacks to swplot figure
    plot            - plots objects to swplot figure
    tooltip         - creates tooltip axis on swplot figure
    transform       - transform objects on swplot figure
    translate       - translate objects on swplot figure
    view            - control 3D view of swplot
    zoom            - zooming objects on swplot figure
 
  Plotting SpinW object on swplot figure:
 
    plotatom        - plots crystal structure
    plotbase        - plots the edges of unit cells on swplot figure
    plotbond        - plots magnetic bonds
    plotcell        - plots the edges of unit cells on swplot figure
    plotchem        - plots polyhedra or chemical bonds
    plotion         - plots magnetic ion properties
    plotmag         - plots magnetic structure
 
  Other related functions:
 
    color           - generates RGB code from color name string
    icomesh         - creates mesh by subdividing icosahedron faces
    logo            - creates the logo of SpinW and saves to a .png file
    patchfacefcn    - callback function for patch face selection
    raytriangle     - finds if a ray crosses a triangle
    setrangegui     - shows a window to change the plotting range of an swplot figure
    subfigure       - changes position of figure window on the screen
    tooltipcallback - call tooltip on clicking an object on an swplot figure
    tooltipstring   - generate tooltip string from the data of a graphical object
 