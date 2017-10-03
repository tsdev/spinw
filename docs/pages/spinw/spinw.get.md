---
{title: spinw.get method, link: spinw.get, summary: GET    Get MATLAB object properties.,
  keywords: sample, sidebar: sw_sidebar, permalink: spinw_get, folder: spinw, mathjax: 'true'}

---

### Syntax

`    handles, get returns an m-by-1 cell array of values, where m is equal`

### Description

    cell array of strings containing property names, GET returns an M-by-N
    cell array of values.  For non-scalar H, if 'PropertyName' is a 
    dynamic  property, GET returns a value only if the property exists in 
    all objects of the array.
  
    V = GET(H, 'InexactPropertyName') returns the value of the specified
    property for the MATLAB object with handle H. GET matches partial and 
    case-insensitive names that are not ambiguous. Inexact name matching 
    applies only to class properties. Dynamic properties require exact name matches.
 
    V = GET(H) returns a structure in which each field name is the name of
    a user-gettable property of H and each field contains the value of that
    property.  If H is non-scalar, GET returns a struct array with 
    dimensions M-by-1, where M = numel(H).  If H is non-scalar, GET does 
    not return dynamic properties.
 
    GET(H) displays the names of all user-gettable properties and their 
    current values for the MATLAB object with handle H.  The class can 
    override the GETDISP method to control how this information is 
    displayed.  H must be scalar.
 

### See Also

[also] \| [get] \| [spinw](spinw) \| [spinw/getdisp](spinw_getdisp) \| [handle]
Help for spinw/get is inherited from superclass MATLAB.MIXIN.SETGET
    Reference page in Doc Center
       doc spinw/get

{% include links.html %}
