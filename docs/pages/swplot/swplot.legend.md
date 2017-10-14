---
{title: swplot.legend, link: swplot.legend, summary: adds legend to the swplot figure,
  keywords: sample, sidebar: sw_sidebar, permalink: swplot_legend, folder: swplot,
  mathjax: 'true'}

---
  
### Syntax
  
`swplot.legend`
  
`swplot.legend(switch, hFigure)`
 
`status = swplot.legend`
 
### Description
  
`swplot.legend` adds legend to the active swplot figure.
   
`swplot.legend(switch, hFigure)` adds/removes/refreshes the legend on the
swplot figure referenced by the `hFigure` handle depending on the
`switch` string.
 
### Examples
  
This example shows how the default legend for arrow and circle objects
looks like.
 
```matlab
swplot.plot('type','arrow','position',rand(3,10,2)*10-5,'legend',1,'color','gold')
swplot.plot('type','circle','position',rand(3,10,2)*10-5,'R',1,'legend',1,'color','purple')
swplot.zoom
swplot.legend
```
 
{% include image.html file="generated/swplot_l_1.png" alt="swplot.legend" %}
  
### Input Arguments
  
`switch`
: One of the following string:
  * `'on'`                show legend,
  * `'off'`               hide legend,
  * `'refresh'`           redraw legend,
  * `'-'`\|`'--'`\|`'none'` change the linestyle of the legend frame.
 
  Default value is `'on'`.
  
`hFigure`
: Handle of the swplot figure, default value is the handle of the active
  figure.
 

{% include links.html %}
