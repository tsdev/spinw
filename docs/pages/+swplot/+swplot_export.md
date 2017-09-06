---
{title: export( ), keywords: sample, summary: exports swplot figure into raster/vector image,
  sidebar: sw_sidebar, permalink: +swplot_export.html, folder: +swplot, mathjax: 'true'}

---
  exports swplot figure into raster/vector image
 
  SWPLOT.EXPORT('Option1',Value1,...)
 
  The function will remove the tooltip before exporting to raster image.
  For vector graphics, also the legend will be removed as it causes a bug
  in the Matlab print command. Also vector graphics export does not support
  transparency, thus all transparency will be removed from the figure. Be
  careful, the vector image filesize can be quite large if there are many
  object on the plot. To reduce the file size, try reducing the npatch and
  nmesh values to reduce the number of faces per object. The function uses
  the Matlab built-in print() command after preparing the figure. All
  figure property restored after export.
 
  Options:
 
  figure    Handle of the swplot figure. Default is the selected figure.
  filename  String, name of the image file. Image type will be determined
            based on the extension. Supported graphics formats:
                .png    Raster image.
                .eps    Vector image.
            If no filename provided, 
  res       Resolution for raster images in dpi, default value is 300. Set
            it to 0, to save the image with the screen resolution.
 
  See also PRINT.
 
