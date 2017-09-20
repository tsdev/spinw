---
{title: spinw.export method, link: spinw.export, summary: export data from spinw object
    into different file formats, keywords: sample, sidebar: sw_sidebar, permalink: spinw_export.html,
  folder: spinw, mathjax: 'true'}

---

### Syntax

`out = export(obj,Name,Value)`

### Description

Different part of the spinw object data can be exported selected by the
'format' option. Right now the following formats are supported:
 
'pcr'     Creates part of a .pcr file used by FullProf. It exports the
          atomic positions.
'spt'     Creates a Jmol script, that reproduce the same plot as used the
          built in spinw.plot() function. Any additional parameter of the
          spinw.plot() function can be used.
'MC'      Exports a proprietary file format for Monte Carlo simulations.
 
 
Other general options:
 
path      Path to a file into which the data will be exported, 'out' will
          be true if the file succesfully saved, otherwise false.
fid       File identifier that is already opened in Matlab using the
          fid = fopen() function. 'out' will be the input fid. Don't
          forget to close the text file afterwards.
 
 
Format related options:
 
PCR format:
perm      Permutation of the xyz atomic positions, default value is [1 2 3].
 
MC format:
boundary  Boundary conditions of the extended unit cell.
                'free'  Free, interactions between extedned unit cells are
                        omitted.
                'per'   Periodic, interactions between extended unit cells
                        are retained.
            Default is {'per' 'per' 'per'}.
 
If neither 'path' nor 'fid' is given, the 'out' will be a cell containing
strings for each line of the text output.
 

### Examples

cryst = sw('test.cif');
cryst.export('format','pcr','path','test.pcr');
In this example the crystal structure is imported from the test.cif file,
and the atomic positions are saved into the test.pcr file for FullProf
refinement (the pcr file needs additional text to work with FullProf).
Links:
Jmol Wiki: http://wiki.jmol.org/index.php/Main_Page
FullProf:  https://www.ill.eu/sites/fullprof

