---
{title: swplot.patchfacefcn, link: swplot.patchfacefcn, summary: callback function
    for patch face selection, keywords: sample, sidebar: sw_sidebar, permalink: swplot_patchfacefcn,
  folder: swplot, mathjax: 'true'}

---
  
### Syntax
  
`patchfacefcn(hPatch, hit, callbackFun, selection)`
  
`patchfacefcn(hPatch, hit, callbackFun, selection, dLim)`
 
### Description
  
`patchfacefcn(hPatch, hit, callbackFun, selection)` finds the index of the
face in a patch object which was clicked on by the mouse. The function
should be used as a callback function for the `ButtonDownFcn` event of
the patch object and it will call a user defined function with the same
arguments as the `ButtonDownFcn` call, plus adding an extra argument, the
face index. Thus the user defined callback function needs to have the
following header:
``` matlab
callbackFun(hPatch,hit,faceIndex)
```
   
The function can detect if the mouse click was on a face or on an edge of
the patch object.
   
`patchfacefcn(hPatch, hit, callbackFun, selection, dLim)` adds optional
control on the selectivity whether the click was on a face of a patch
object.
 
### Examples
  
The color of any face of the red triangulated icosahedron will be
changed from red to green if clicked on with the mouse.
 
```matlab
mesh = swplot.icomesh(1)
V = mesh.X
F = mesh.Triangulation
hPatch = patch('Faces',F,'Vertices',V,'FaceColor','r','EdgeColor','none')
axis equal
axis off
box on
view(3)
camlight('right')
hPatch.FaceColor = 'flat'
hPatch.FaceVertexCData = repmat([1 0 0],[size(F,1) 1])
fun = @(hPatch,hif,face)set(hPatch,'FaceVertexCData',[hPatch.FaceVertexCData(1:(face-1),:); [0 1 0]; hPatch.FaceVertexCData((face+1):end,:)])
hPatch.ButtonDownFcn = @(hPatch,hit)swplot.patchfacefcn(hPatch,hit,fun,'face')
```
  
### Input Arguments
  
`hPatch`
: Handle of the patch object.
  
`hit`
: Hit object that contains the point where the object was hit.
  
`callbackFun`
: User function that will be called in case of a click event
    on `hPatch` object. It should have the following header:
        `callbackFun(hPatch,hit,faceIndex)`
    where `faceIndex` contains the index of the face that was
    clicked on, it can contain a single index or more depending
    on the selection type.
  
`selection`
: String, that defines three diferent selection criteria when the
  `callbackfun()` function will be called:
  * `'face'`  The `callbackFun()` will be triggered if a face
              was clicked (`numel(faceIndex)==1`).
  * `'edge'`  The `callbackfun()` will be triggered if an edge
              is clicked (`numel(faceIndex)==2`).
  * `'all'`   The `callbackfun()` will be triggered for both
              faces and edges.
  
`dLim`
: Upper limit of the absolute value of the determinant that
  determines whether a point is on a plane spanned by the two
  edges of a triangle. Default value is $$10^{-7}$$ (tested).
  
### See Also
  
[patch](https://www.mathworks.com/help/matlab/ref/patch.html)
 

{% include links.html %}
