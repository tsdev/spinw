---
title: getpref( )
keywords: sample
summary: "returns SpinW global preferences"
sidebar: product1_sidebar
permalink: getpref.html
folder: +swpref
mathjax: true
---
  returns SpinW global preferences
 
  rPref = swpref.getpref
 
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
 
  See also SWPREF.SETPREF.
