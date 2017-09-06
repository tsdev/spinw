---
{title: swprefsetpref( ), keywords: sample, summary: sets SpinW global preferences,
  sidebar: sw_sidebar, permalink: +swpref_setpref.html, folder: +swpref, mathjax: 'true'}

---
sets SpinW global preferences
 
swpref.setpref(prefName, value)
 
The preferences are reset after every restart of Matlab, unlike the
Matlab built-in preferences that are persistent between Matlab sessions.
If you want certain preferences to keep after closing matlab, define them
in the <a href="matlab:edit('startup.m')">startup.m</a> file.
 
swpref.setpref() sets the value of the prefName in the SpinW global
preferences.
 
swpref.setpref('default')
 
Resets all preference values to the default one.
 
See also SWPREF.SETPREF.
 
