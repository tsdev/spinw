---
{title: swplot.mouse, link: swplot.mouse, summary: adds mouse callbacks to swplot
    figure, keywords: sample, sidebar: sw_sidebar, permalink: swplot_mouse.html, folder: swplot,
  mathjax: 'true'}

---
 
SWPLOT.MOUSE({hFigure},{perspective})
 
Adds rotation and zoom functionality to any swplot figure. The following
mouse actions are supported:
  - mouse-drag        Rotation of objects.
  - ctrl+mouse-drag   Shift of objects (pan).
  - mouse-wheel       Zoom of objects.
  - ctrl+mouse-wheel  Change perspective and switch to perspective
                      projection.
 
Input:
 
hFigure       Handle of the swplot figure. Default is the selected
              figure.
perspective   String determines whether camera projection mode is changed
              automatically between orthographic (zooming withouth CTRL 
              key pressed) and perspective (zooming with CTRL key
              pressed):
                  'auto'      Automatic switching (default).
                  'fix'       No switching.
 
See also CAMPROJ.
 

