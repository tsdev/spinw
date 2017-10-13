---
{title: swplot.getdata, link: swplot.getdata, summary: gets the data stored in an
    swplot figure, keywords: sample, sidebar: sw_sidebar, permalink: swplot_getdata,
  folder: swplot, mathjax: 'true'}

---
  
### Syntax
  
`data = swplot.getdata`
  
`data = swplot.getdata(hFigure)`
 
`data = swplot.getdata(field)`
 
### Description
  
`data = swplot.getdata` gets all the object data stored in the active
swplot figure.
 
`data = swplot.getdata(hFigure)` get all object data stored in the swplot
figure identified by the `hFigure` handle.
 
`data = swplot.getdata(field)` loads only the given field of the data
structure.
  
### Examples
 
This example shows how the data of all objects on a 3D SpinW plot can be
retrieved.
 
```matlab
model = sw_model('triAF',1)
plot(model)
swplot.getdata
```
*Output*
```
  struct with fields:
          h: [1×1 Transform]
     button: [1×1 struct]
       axis: [1×1 Axes]
     legend: [1×1 struct]
      light: [1×1 Light]
       icon: [1×1 struct]
       base: [3×3 double]
    tooltip: [1×1 struct]
    objects: [1×44 struct]
        obj: [1×1 spinw]
      range: [1×1 struct]
      param: {}
```
 
  
### Input Arguments
  
`hFigure`
: Handle of the swplot figure, default value is the active figure.
  
`field`
: String, determines the requested field name. If omitted, all
  stored fields are returned.
  
### See Also
  
[getappdata](https://www.mathworks.com/help/matlab/ref/getappdata.html)
 

{% include links.html %}
