---
title: pref( )
keywords: sample
summary: "returns SpinW global preferences"
sidebar: product1_sidebar
permalink: pref.html
folder: +swpref
mathjax: true
---
  returns SpinW global preferences
 
  rPref = swpref.getpref
 
  The preferences are reset after every restart of Matlab, unlike the
  Matlab built-in preferences that are persistent between Matlab sessions.
  If you want certain preferences to keep after closing matlab, define them
  in the <a href="matlab:edit('startup.m')">startup.m</a> file.
 
  swpref.getpref() returns the names, values and labels of each
  preferences. Default values are returned, if no values are saved. rPref
  is a struct with field names 'name', 'label' and 'val'. Each field is a
  cell.
 
  rPref = swpref.getpref(pName, {simple})
 
  Returns only the requested SpinW preference name, value and label in a
  struct. Each field contains the requested value. If a second argument is
  given (simple) with any value, only the value of the preference is
  returned.
 
  rPref = swpref.getpref('default')
 
  Returns the default names, values and labels of each preferences.
 
  See also GETPREF, SETPREF, SWPREF.SETPREF.
