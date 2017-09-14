---
{title: spinw.getdisp method, link: spinw.getdisp, summary: GETDISP    Specialized
    MATLAB object property display., keywords: sample, sidebar: sw_sidebar, permalink: spinw_getdisp.html,
  folder: spinw, mathjax: 'true'}

---
   GETDISP is called by GET when GET is called with no output argument 
   and a single input parameter H an array of handles to MATLAB objects.  
   This method is designed to be overridden in situations where a
   special display format is desired to display the results returned by
   GET(H).  If not overridden, the default display format for the class
   is used.
 
   See also spinw, spinw/GET, handle
Help for spinw/getdisp is inherited from superclass MATLAB.MIXIN.SETGET
   Reference page in Doc Center
      doc spinw/getdisp

