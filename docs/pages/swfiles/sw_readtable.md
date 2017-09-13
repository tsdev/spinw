---
{title: sw_readtable( ), summary: reads tabular data, keywords: sample, sidebar: sw_sidebar,
  permalink: sw_readtable.html, folder: swfiles, mathjax: 'true'}

---
reads tabular data
 
dat = SW_READTABLE(dataSource, {delimiter},{nHeader})
 
Function reads tabular data that has arbitrary header lines denoted with
# and the last header line is followed by a column name line. The data
can be arbitrary combination of strings and numbers. A predefined field
will be added to the imported data 'MODE', this field contains any
additional string given between data rows, otherwise empty.
 
Input:
 
dataSource    Data source, can be file, url or string.
delimiter     Delimiter of the data, default is whitespace.
 
Output:
 
dat       Structure with field defined in the data.
 
Example:
 
The following file (test.dat) is given as input:
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
 
The command to import the data
>> dat = sw_readtable('test.dat');
 
The dat variable will contain the fields MODE, Q, ENlim, I, EN and s and
it will have 7 entry. The mode variable will contain '[Mxx] [1 0 0]'
string for the first 4 entry and '[Myy] [1 0 0]' for the last 3 entry.
For example the field Q has 3 elements per entry, to extract all Q points
into a matrix use the command:
>> Q = reshape([dat(:).Q],3,[])';
 
