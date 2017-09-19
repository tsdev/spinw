---
{title: sw_update( ), link: sw_update, summary: updates the SpinW installation from
    the internet, keywords: sample, sidebar: sw_sidebar, permalink: sw_update.html,
  folder: swfiles, mathjax: 'true'}

---

### Syntax

` `

### Description

 
sw_update creates a new folder with the latest release beside the current
SpinW installation and add the new version to the working path (and
removing the old one).
 

### Input Arguments

% `installDir`
:r    Folder name, where the new version is installed. Default is

% `the`
:t folder of the current version of SpinW. If

% `installDir`
:r == '.' update will be installed to current

% ``
:

### Output Arguments

onlineVer     If output is expected, the revision number of the online
SpinW is given. Optional.

