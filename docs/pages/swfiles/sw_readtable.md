---
{title: sw_readtable, link: sw_readtable, summary: reads tabular data from text, keywords: sample,
  sidebar: sw_sidebar, permalink: sw_readtable, folder: swfiles, mathjax: 'true'}

---
  
### Syntax
  
`dat = sw_readtable(datasource)`
  
`dat = sw_readtable(datasource,delimiter,nheader)`
 
### Description
  
`dat = sw_readtable(datasource)` reads tabular data and converts it into
a struct with as many number of elements as many data rows can be found
in the source. The output field names are determined by the header line
that begins with a string.  Moreover multiple data columns can be
combined using the same column name in the header line multiple times and
an additional index in brackets.
 
The reader also supports tagging of row of data with a string. A line
that begins with `#` defines a tag which will be applied to all following
data until a new tag is defined. The tags are saved in the `tag` field
of the output structure.
  
### Examples
  
The content of the input file (`test.dat`) is the following:
 
```
# TEST DATA
Q(1) Q(2)        Q(3) ENlim(1) ENlim(2) I(1)  EN(1)  s(1) I(2)   EN(2)   s(2)
# [Mxx] [1 0 0]
0     1        2.9992   0       15      1    3.7128   1.0   1   8.6778    1.0
0     1        2.8993   0       15      1    7.0000   1.0   1   11.1249   1.0
0     1        2.7993   0       20      1   13.8576   1.0   0   0.0       0.0
0     1        2.6994   0       20      1   17.3861   1.0   0   0.0       0.0
# [Myy] [1 0 0]
0     1.0000   2.0000   0       25      1   20.2183   1.0   0   0.0       0.0
0     1.1000   2.0000   15      30      1   22.7032   1.0   0   0.0       0.0
0     1.2000   2.0000   20      35      1   25.1516   1.0   0   0.0       0.0
```
 
The command to import the data:
 
```matlab
dat = sw_readtable('test.dat')
dat(1)
```
*Output*
```
  struct with fields:
      tag: '# [Mxx] [1 0 0]'
        Q: [0 1 2.9992]
    ENlim: [0 15]
        I: [1 1]
       EN: [3.7128 8.6778]
        s: [1 1]
```
 
```matlab
Q = reshape([dat(:).Q],3,[])'
```
*Output*
```
Q =
         0    1.0000    2.9992
         0    1.0000    2.8993
         0    1.0000    2.7993
         0    1.0000    2.6994
         0    1.0000    2.0000
         0    1.1000    2.0000
         0    1.2000    2.0000
```
 
 
Here the imported `dat` variable will contain the fields `tag`, `Q`,
`ENlim`, `I`, `EN` and `s` and it will have 7 entry. The `tag` variable
will contain the `'[Mxx] [1 0 0]'` string for the first 4 entries and
`'[Myy] [1 0 0]'` for the last 3 entries. For example the field `Q` has 3
elements per entry and the last command above extracts all $$Q$$ points
into a matrix.
  
### Input Arguments
  
`dataSource`
: Data source, can be file, url or string (must contain the newline
  character).
  
`delimiter`
: Delimiter of the data, default value is whitespace.
 
`nheader`
: Number of header lines to be skipped in the beginning of the file.
  Default value is 0.
  
### Output Arguments
  
`dat`
: Struct with fields defined in the data.
 

{% include links.html %}
