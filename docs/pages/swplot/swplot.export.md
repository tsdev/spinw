---
{title: swplot.export, link: swplot.export, summary: exports swplot figure into raster/vector
    image, keywords: sample, sidebar: sw_sidebar, permalink: swplot_export, folder: swplot,
  mathjax: 'true'}

---
  
### Syntax
  
`swplot.export(Name,Value)`
  
### Description
  
`swplot.export(Name,Value)` exports an swplot figure into a raster/vector
image. The function will remove the tooltip before exporting to raster
image. For vector graphics, also the legend will be removed as it causes
a bug in the Matlab [print](https://www.mathworks.com/help/matlab/ref/print.html) command. Also, vector graphics export
does not support transparency, thus all transparency will be removed from
the figure. Be careful, the vector image filesize can be quite large if
there are many object on the figure. To reduce the file size, try
reducing the $$n_{patch}$$ and $$n_{mesh}$$ values to reduce the number of
faces per object. The function uses the Matlab built-in [print](https://www.mathworks.com/help/matlab/ref/print.html)
command after preparing the figure. All figure property restored after
export.  
  
### Name-Value Pair Arguments
  
`'figure'`
: Handle of the swplot figure. Default value is the active figure.
  
`'filename'`
: String, name of the image file. Image type will be determined
  based on the extension. Supported graphics formats:
  * `.png`    Raster image.
  * `.eps`    Vector image.
 
  If no filename provided, the function returns without printing.
  
`'res'`
: Resolution for raster images in dpi, default value is 300. Set
  it to 0, to save the image with the screen resolution.
  
### See Also
  
[print](https://www.mathworks.com/help/matlab/ref/print.html)
 

{% include links.html %}
