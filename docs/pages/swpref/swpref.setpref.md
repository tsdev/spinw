---
{title: swpref.setpref, link: swpref.setpref, summary: sets SpinW global preferences,
  keywords: sample, sidebar: sw_sidebar, permalink: swpref_setpref, folder: swpref,
  mathjax: 'true'}

---
  
### Syntax
  
`swpref.setpref(prefName, value)`
  
`swpref.setpref('default')`
 
### Description
  
`swpref.setpref(prefName, value)` sets the value of `prefName`
preferences.
 
`swpref.setpref('default')` resets all preference values to the default one.
 
{% include note.html content=" The preferences are reset after every restart of Matlab, unlike the
Matlab built-in preferences that are persistent between Matlab sessions.
If you want certain preferences to keep after closing matlab, define them
in the `startup.m` file." %}
  
### See Also
  
[swpref.getpref](swpref_getpref)
 

{% include links.html %}
