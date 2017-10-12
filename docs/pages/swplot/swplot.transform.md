---
{title: swplot.transform, link: swplot.transform, summary: transform objects on swplot
    figure, keywords: sample, sidebar: sw_sidebar, permalink: swplot_transform, folder: swplot,
  mathjax: 'true'}

---

### Syntax

`swplot.transform(m, {hfigure})`

### Description

Transforms the objects on the active swplot figure using the
transformation matrix M.
 

### Input Arguments

`M`
:   Transformation matrix with possible dimensions:
        4x4     This follows the Matlab standard for hgtransform.
        3x4     This is the SpinW format for space group 
                transformations. 
        3x3     This defines the rotation matrix only.
    Setting M to 0 returns to the plot to the original orientation
    (equivalent to M=eye(4)).

`hFigure`
:   Handle of the swplot figure window, optional.

`M`
:LOT.TRANSFORM({hFigure})

`Returns`
: the transformation matrix (4x4) of the active swplot figure.

`If`
:figure is created without the hgtransform object, the

`transformation`
:rmation matrix moves the camera.

### See Also

[swplot.figure](swplot_figure) \| [hgtransform]

{% include links.html %}
