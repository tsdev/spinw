---
{title: swplot.findobj, link: swplot.findobj, summary: finds object data on swplot
    figure, keywords: sample, sidebar: sw_sidebar, permalink: swplot_findobj, folder: swplot,
  mathjax: 'true'}

---
  
### Syntax
  
`sObj = swplot.findobj(Name,Value)`
 
`sObj = swplot.findobj(hFigure,Name,Value)`
  
### Description
  
`sObj = swplot.findobj(Name,Value)` finds graphical objects on the active
swplot figure hFigure which have the given property name-value pairs. The
possible property names are:
* `handle`    Handle of the graphical object.
* `objID`     Unique number of the object (increasing integer numbers).
* `name`      Name of the object, identifies groups, such as `'atom'` for
              all atoms.
* `label`     Label of the objects, can identify types of atoms, etc.
              it will accept sub strings, e.g. `'Cr'` parameter would
              match both `'Cr1 Cr3+'` and `'Cr2 Cr3+'` labels.
 
`sObj = swplot.findobj(hFigure,Name,Value)` search for objects on the
swplot figure identified by the `hFigure` handle.
 
### Output Arguments
 
`sObj`
: Struct that contains all the data of the found objects.
   
### See Also
 
[swplot.delete](swplot_delete)
 

{% include links.html %}
