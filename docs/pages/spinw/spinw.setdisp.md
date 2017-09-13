---
{title: spinw.setdisp( ), summary: SETDISP    Specialized MATLAB object property display.,
  keywords: sample, sidebar: sw_sidebar, permalink: spinw_setdisp.html, folder: spinw,
  mathjax: 'true'}

---
SETDISP    Specialized MATLAB object property display.
   SETDISP is called by SET when SET is called with no output argument 
   and a single input parameter H an array of handles to MATLAB objects.  
   This method is designed to be overridden in situations where a
   special display format is desired to display the results returned by
   SET(H).  If not overridden, the default display format for the class
   is used.
 
   See also setdisp, [spinw](spinw.html), [spinw/set](spinw_set.html), handle
Help for [spinw](spinw.html)/setdisp is inherited from superclass MATLAB.MIXIN.SETGET
   Reference page in Doc Center
      doc [spinw](spinw.html)/setdisp

