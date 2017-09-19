---
{title: swplot.patchfacefcn, link: swplot.patchfacefcn, summary: callback function
    for patch face selection, keywords: sample, sidebar: sw_sidebar, permalink: swplot_patchfacefcn.html,
  folder: swplot, mathjax: 'true'}

---

### Syntax

` `

### Description

clicked on by the mouse. The function should be used as a callback
function for the ButtonDownFcn event of patch object and it will call a
user defined function with the same arguments as the ButtonDownFcn call,
but adding an extra argument, the face index. Thus the user defined
callback function will have the following header:
 
callbackfun(obj,hit,faceIndex)
 
The function can detect if the mouse click was on a face or on an edge of
the patch object.
 

### Examples



### Input Arguments

% `obj`
: object.

% `hit`
:   Hit object that contains the point where the object was hit.

% `callbackfun`
:   User function that will be called in case of a click event
  It should have the following header:
 lbackfun(obj,hit,faceIndex)
 ace Index contains the index of the face that was
  on, it can contain a single index or more depending
 selection type.

% `selection`
:   String defines three diferent selection criteria when the
 kfun() function will be called:
 ce'  The callbackfun() will be triggered if a face
      was clicked (numel(faceIndex)==1).
 ge'  The callbackfun() will be triggered if an edge
      is clicked (numel(faceIndex)==2).
 l'   The callbackfun() will be triggered for both
      faces and edges.

% `{dLim}`
:   Upper limit of the absolute value of the determinant that
 nes whether a point is on a plane spanned by the two
 f a triangle. Default value is 1e-7 (tested).

### See Also

[patch]

