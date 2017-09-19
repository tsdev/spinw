---
{title: swplot.mouse, link: swplot.mouse, summary: adds mouse callbacks to swplot
    figure, keywords: sample, sidebar: sw_sidebar, permalink: swplot_mouse.html, folder: swplot,
  mathjax: 'true'}

---

### Syntax

` `

### Description

mouse actions are supported:
- mouse-drag        Rotation of objects.
- ctrl+mouse-drag   Shift of objects (pan).
- mouse-wheel       Zoom of objects.
- ctrl+mouse-wheel  Change perspective and switch to perspective
                    projection.
 

### Input Arguments

% `hFigure`
: Handle of the swplot figure. Default is the selected
 igure.

% `perspective`
: String determines whether camera projection mode is changed
 utomatically between orthographic (zooming withouth CTRL 
 ey pressed) and perspective (zooming with CTRL key
 ressed):
    'auto'      Automatic switching (default).
    'fix'       No switching.

### See Also

[camproj]

