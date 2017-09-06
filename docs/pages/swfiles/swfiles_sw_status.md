---
{title: sw_status( ), keywords: sample, summary: timer function that displays also the remaining time,
  sidebar: sw_sidebar, permalink: swfiles_sw_status.html, folder: swfiles, mathjax: 'true'}

---
  timer function that displays also the remaining time
 
  SW_STATUS(percent, {mode},{fid},{title})
 
  Input:
 
  percent   Percentage of the calculation that is done.
  mode      Determines the time estimation, optional parameter:
                1   Starts the time estimation.
                0   Displays of the remaining time. (default)
                2   Calculation finished.
  fid       File identifier to print the output:
                0   Do nothing.
                1   Text output to the Command Window. Default.
                2   Graphical output, using the waitbar() function.
 
  See also WAITBAR.
 
