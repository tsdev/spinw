---
{title: sw_update( ), summary: updates the SpinW installation from the internet, keywords: sample,
  sidebar: sw_sidebar, permalink: sw_update.html, folder: swfiles, mathjax: 'true'}

---
updates the SpinW installation from the internet
 
SW_UPDATE()
 
{onlineRev} = SW_UPDATE(installDir)
 
sw_update creates a new folder with the latest release beside the current
SpinW installation and add the new version to the working path (and
removing the old one).
 
Input:
 
installDir    Folder name, where the new version is installed. Default is
              the parent folder of the current version of SpinW. If
              installDir == '.' update will be installed to current
              folder.
 
Output:
 
onlineVer     If output is expected, the revision number of the online
              SpinW is given. Optional.
 
