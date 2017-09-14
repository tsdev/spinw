---
{title: spinw.set( ), summary: SET    Set MATLAB object property values., keywords: sample,
  sidebar: sw_sidebar, permalink: spinw_set.html, folder: spinw, mathjax: 'true'}

---
   SET(H,'PropertyName',PropertyValue) sets the value of the specified 
   property for the MATLAB object with handle H.  If H is an array of 
   handles, the specified property's value is set for all objects in H.  
 
   SET(H,'InexactPropertyName',PropertyValue) sets the value of the specified 
   property for the MATLAB object with handle H. SET matches partial and 
   case-insensitive names that are not ambiguous. Inexact name matching 
   applies only to class properties. Dynamic properties require exact name matches.
 
   SET(H,'PropertyName1',Value1,'PropertyName2',Value2,...) sets multiple
   property values with a single statement. 
 
   SET(H,pn,pv) sets the named properties specified in the cell array of
   strings pn to the corresponding values in the cell array pv for all
   objects specified in H.  The cell array pn must be 1-by-N, but the cell
   array pv can be M-by-N where M is equal to length(H), so that each
   object will be updated with a different set of values for the list of
   property names contained in pn.
 
   Given S a structure whose field names are object property names, 
   SET(H,S) sets the properties identified by each field name of S with
   the values contained in the structure.
 
   Note that it is permissible to use property/value string pairs, 
   structures, and property/value cell array pairs in the same call to
   SET.
 
   A = SET(H, 'PropertyName') returns the possible values for the 
   specified property of the object with handle H.  The returned array
   is a cell array of possible value strings or an empty cell array if
   the property does not have a finite set of possible string values.
 
   SET(H,'PropertyName') displays the possible values for the specified
   property of object with handle H.
 
   A = SET(H) returns the names of the user-settable properties and their
   possible values for the object with handle H.  H must be scalar.  
   The return value is a  structure whose field names are the names of the
   user-settable properties of H, and whose values are cell arrays of
   possible property values or empty cell arrays.
 
   SET(H) displays the names and possible values for all user-settable
   properties of scalar object H.  The class can override the SETDISP 
   method to control how this information is displayed. 
 
   See also SET, spinw, spinw/SETDISP, HANDLE
Help for spinw/set is inherited from superclass MATLAB.MIXIN.SETGET
   Reference page in Doc Center
      doc spinw/set

