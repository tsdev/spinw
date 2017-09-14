---
{title: swplot.patchfacefcn( ), summary: callback function for patch face selection,
  keywords: sample, sidebar: sw_sidebar, permalink: swplot_patchfacefcn.html, folder: swplot,
  mathjax: 'true'}

---
 
PATCHFACEFCN(obj,hit,callbackfun,selection, {dLim})
 
The function can find the index of the face in a patch object which was
clicked on by the mouse. The function should be used as a callback
function for the ButtonDownFcn event of patch object and it will call a
user defined function with the same arguments as the ButtonDownFcn call,
but adding an extra argument, the face index. Thus the user defined
callback function will have the following header:
 
          callbackfun(obj,hit,faceIndex)
 
The function can detect if the mouse click was on a face or on an edge of
the patch object.
 
Input:
 
obj           The patch object that calls the function.
hit           Hit object that contains the point where the object was hit.
callbackfun   User function that will be called in case of a click event
              on obj. It should have the following header:
                  callbackfun(obj,hit,faceIndex)
              where face Index contains the index of the face that was
              clicked on, it can contain a single index or more depending
              on the selection type.
selection     String defines three diferent selection criteria when the
              callbackfun() function will be called:
                  'face'  The callbackfun() will be triggered if a face
                          was clicked (numel(faceIndex)==1).
                  'edge'  The callbackfun() will be triggered if an edge
                          is clicked (numel(faceIndex)==2).
                  'all'   The callbackfun() will be triggered for both
                          faces and edges.
{dLim}        Upper limit of the absolute value of the determinant that
              determines whether a point is on a plane spanned by the two
              edges of a triangle. Default value is 1e-7 (tested).
 
Example:
 
See also PATCH.
 

