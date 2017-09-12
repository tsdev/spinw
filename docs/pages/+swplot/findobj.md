---
{title: swplot.findobj( ), summary: finds object data on swplot figure, keywords: sample,
  sidebar: sw_sidebar, permalink: swplot_findobj.html, folder: +swplot, mathjax: 'true'}

---
finds object data on swplot figure
 
sObj = swplot.findobj(hFigure,'PropertyName',value,...)
 
Finds objects on hFigure that all have the given property values. The
possible property names are:
 
  handle      Handle of the graphical object.
  number      Unique number of the object (increasing integer numbers).
  name        Name of the object, identifies groups, such as 'atom' for
              all atoms.
  label       Label of the objects, can identify types of atoms, etc. if
              will accept sub strings. E.g. 'Cr' option would match 'Cr1
              Cr3+' and 'Cr2 Cr3+' labels.
 
sObj is a struct that contains all the stored data corresponding to the
found objects.
 
sObj = swplot.findobj('PropertyName',value,...)
 
Finds object on the active swplot figure.
 
