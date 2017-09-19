---
{title: swplot.transform, link: swplot.transform, summary: transform objects on swplot
    figure, keywords: sample, sidebar: sw_sidebar, permalink: swplot_transform.html,
  folder: swplot, mathjax: 'true'}

---

### Syntax

` `

### Description

transformation matrix M.
 

### Input Arguments

% `M`
:  Transformation matrix with possible dimensions:
     This follows the Matlab standard for hgtransform.
     This is the SpinW format for space group 
     transformations. 
     This defines the rotation matrix only.

% `Setting`
:M to 0 returns to the plot to the original orientation

% `(equivalent`
:ent to M=eye(4)).

% `hFigure`
:  Handle of the swplot figure window, optional.

% `M`
:OT.TRANSFORM({hFigure})

% `Returns`
:the transformation matrix (4x4) of the active swplot figure.

% `If`
:igure is created without the hgtransform object, the

% `transformation`
:mation matrix moves the camera.

### See Also

[swplot.figure](swplot_figure.html) and [hgtransform]

