---
{title: disp( ), keywords: sample, summary: prints the spinw data structure in readable format onto the Command Window,
  sidebar: sw_sidebar, permalink: '@spinw_disp.html', folder: '@spinw', mathjax: 'true'}

---
  prints the spinw data structure in readable format onto the Command Window
 
  {swDescr} = DISPLAY(obj)
 
  Input:
 
  obj       spinw class object.
 
  Output:
 
  swdescr   If output variable is given, the description of the obj object
            will be output into the swdescr variable, instead of being
            written onto the Command Window/file. Optional.
 
  Example:
 
  crystal = spinw;
  swFields = display(crystal);
 
  See also SPINW.
 
