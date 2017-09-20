---
{title: swpref.getpref, link: swpref.getpref, summary: returns SpinW global preferences,
  keywords: sample, sidebar: sw_sidebar, permalink: swpref_getpref.html, folder: swpref,
  mathjax: 'true'}

---

### Syntax

`rpref = swpref.getpref`

### Description

The preferences are reset after every restart of Matlab, unlike the
Matlab built-in preferences that are persistent between Matlab sessions.
If you want certain preferences to keep after closing matlab, define them
in the <a href="matlab:edit('startup.m')">startup.m</a> file.
 
swpref.getpref() returns the names and current values for all
preferences. rPref is a struct.
 
rPref = swpref.getpref(pName, {true})
 
Returns only the requested SpinW preference name, value and label in a
struct. Each field contains the requested value. If a second argument is
given (true), only the value of the preference is returned.
 
rPref = swpref.getpref('default')
 
Returns the default names, values and labels of each preferences.
 

### See Also

[swpref.setpref](swpref_setpref.html)

