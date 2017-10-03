---
{title: sw_status, link: sw_status, summary: timer and remaining time estimator, keywords: sample,
  sidebar: sw_sidebar, permalink: sw_status, folder: swfiles, mathjax: 'true'}

---
  
### Syntax
  
`sw_status(percent, {mode},{tid},{title})`
  
### Description
  
`sw_status(percent, {mode},{fid},{title})` can display remaining time of
a calculation that is run for a fixed number of iterations. It can output
the status both in the Command Window and in a pup up window using
[waitbar].
  
### Input Arguments
  
`percent`
: Percentage of the calculation that is already done.
  
`mode`
: Controls the time estimation, optional parameter:
  * `1` Starts the time estimation.
  * `0` Displays of the remaining time, default value.
  * `2` The calculation finished, show a summary.
  
`tid`
: Determines if the elapsed and required time for the calculation is
  displayed. The default value is determined by the `tid` preference
  stored in [swpref]. The following values are allowed:
  * `0` No timing is displayed.
  * `1` Display the timing in the Command Window.
  * `2` Show the timing in a separat pup-up window.
 
`title`
: The text to show in the pup up window.
  
### See Also
  
[waitbar]

{% include links.html %}
