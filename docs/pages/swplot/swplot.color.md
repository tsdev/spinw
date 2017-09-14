---
{title: swplot.color, link: swplot.color, summary: generates RGB code from color name
    string, keywords: sample, sidebar: sw_sidebar, permalink: swplot_color.html, folder: swplot,
  mathjax: 'true'}

---
 
RGB = SWPLOT.COLOR(cName,{index})
 
Input:
 
cName     String that contains the name of the color, either a single
          character (see <a href="matlab: doc ColorSpec">ColorSpec</a>) or use any HTML color name,
          (see http://www.w3schools.com/html/html_colornames.asp).
          For multiple colors, use a cell containing the strings. The
          name of the colors are stored in the <a href="matlab: edit color.dat">color.dat</a> file.
index     If true, the index of the color in the color.dat file is read.
          Index 1 corresponds to the 9th entry (the first 8 entry has
          already names in Matlab). Default is false.
 
Output:
 
RGB       RGB color code, dimensions are [3 nColor], where
          every value is between 0 and 255.
 
Example:
 
  RGB = SWPLOT.COLOR('LightGray')
 
  the output RGB will be [211; 211; 211].
 

