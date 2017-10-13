---
{title: swplot.delete, link: swplot.delete, summary: deletes objects and corresponding
    data from swplot figure, keywords: sample, sidebar: sw_sidebar, permalink: swplot_delete,
  folder: swplot, mathjax: 'true'}

---
  
### Syntax
  
`swplot.delete(objID)`
 
`swplot.delete(hFigure,objID)`
  
### Description
  
`swplot.delete(objID)` deletes objects and their data that corresponds to
the given unique `objID` (integer number) on the active swplot figure.
 
`swplot.delete(hFigure,objID)` deletes objects from the swplot figure
corresponding to `hFigure` handle.
   
If `objID` equals to 0, all objects will be deleted from the swplot
figure.
   
### See Also
 
[swplot.figure](swplot_figure) \| [swplot.add](swplot_add)
 

{% include links.html %}
