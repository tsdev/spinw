---
{title: sw_filelist, link: sw_filelist, summary: lists spinw objects in the Matlab
    workspace or in a .mat file, keywords: sample, sidebar: sw_sidebar, permalink: sw_filelist,
  folder: swfiles, mathjax: 'true'}

---
  
### Syntax
  
`list = sw_filelist(Name,Value)`
  
### Description
  
 `list = sw_filelist(Name,Value)` lists SpinW objects and calculated
 spectra in the Matlab workspace of in a given .mat file.
  
### Examples
 
After calculating a few spectra, we list them:
 
```matlab
tri = sw_model('triAF',1)
sq  = sw_model('squareAF',1)
specSq  = sq.spinwave({[0 0 0] [1 1 0] 21})
specTri = tri.spinwave({[0 0 0] [1 0 0] 21})
sw_filelist
```
*Output*
```
SpinW variables in the Matlab base workspace:
        VarName                                         Title             Creation Date           Completion Date
            ans                                             -                         -                         -
```
 
 
### Name-Value Pair Arguments
  
`'fName'`
: File path in a string pointing to a .mat file, default value is empty
  when no files are checked.
  
`'sort'`
: Selects the column to sort the generated table with positive/nagetive
  sign means ascending, descending:
  * `'+1'`/`'-1'`    variable name,
  * `'+2'`/`'-2'`    title,
  * `'+3'`/`'-3'`    creation date,
  * `'+4'`/`'-4'`    completion date, default value.
  
### Output Arguments
  
`list`
: Cell of strings, lists each simulation data in the Matlab
  memory, even data stored in cells.
  
### See Also
  
[spinw.anneal](spinw_anneal) \| [spinw.spinwave](spinw_spinwave)
 

{% include links.html %}
