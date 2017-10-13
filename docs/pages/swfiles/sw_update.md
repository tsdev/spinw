---
{title: sw_update, link: sw_update, summary: updates the SpinW installation from the
    internet, keywords: sample, sidebar: sw_sidebar, permalink: sw_update, folder: swfiles,
  mathjax: 'true'}

---
  
### Syntax
  
`sw_update`
  
`onlineRev = sw_update(installDir)`
 
### Description
  
`sw_update` creates a new folder with the latest release beside the
current SpinW installation, downloads the newses SpinW version and adds
the new version to the Matlab search path and also removes the old
version from the search path. If the search path is defined in the
`startup.m` file, it has to be changed manually.
   
Each step of the update process can be controlled by the user via
the interactive Command Line provided by the function.
 
### Input Arguments
  
`installDir`
: Folder name, where the new version is installed. Default is
  the parent folder of the current version of SpinW. If
  `installDir` is `'.'`, the update will be installed to current
  folder.
  
### Output Arguments
  
`onlineVer`  If output is defined, the revision number of the online
              SpinW is given, optional.
 

{% include links.html %}
