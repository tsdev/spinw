---
{title: swplot.color, link: swplot.color, summary: generates RGB code from color name,
  keywords: sample, sidebar: sw_sidebar, permalink: swplot_color, folder: swplot,
  mathjax: 'true'}

---
  
### Syntax
  
`RGB = swplot.color(cName)`
 
`RGB = swplot.color(cName,index)`
  
### Description
  
`RGB = swplot.color(cName)` reads the color RGB values from the
`color.dat` file corresponding to the given color name `cName`. The
color name can be either a single character (see [colorspec](https://www.mathworks.com/help/matlab/ref/colorspec.html)) or
any [HTML color name](https://www.w3schools.com/colors/colors_names.asp)
that is stored in the `color.dat` file.
  
`RGB = swplot.color(cName,index)` if `index` is true, RGB code
corresponding to the `cName` color index is read.
 
### Examples
  
Read the RGB code corresponding to light gray:
```matlab
RGB = swplot.color('LightGray')
```
*Output*
```
RGB =
   211
   211
   211
```
 
  
### Input Arguments
  
`cName`
: String of a color name. For multiple colors, use a cell of strings.
  
`index`
: If `true`, instead of the color name, `cName` means the index of the
  color in the `color.dat` file. index 1 corresponds to the 9th entry
  (the first 8 entry are standard Matlab color names), default value is
  `false`.
  
### Output Arguments
  
`RGB`
: RGB color codes in a matrix with dimensions of $$[3\times n_{color}]$$, where
  every value is an integer between 0 and 255.
 

{% include links.html %}
