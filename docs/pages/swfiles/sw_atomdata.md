---
{title: sw_atomdata( ), link: sw_atomdata, summary: returns information on elements
    stored in the atom.dat file, keywords: sample, sidebar: sw_sidebar, permalink: sw_atomdata.html,
  folder: swfiles, mathjax: 'true'}

---

### Syntax

` `

### Description



### Examples

sw_atomdata('H','radius') = 0.37
If the atom label does not exists, the function returns radius = 1,
color = [255 167 0].
optional second output is 'atomLabel' that contains the name of the atom
clean.

### Input Arguments

% `atomSymb`
:  String of the name of the atom, for example 'He' or the atomic

% `number`
:. If the string contains whitespace character, the

% `second`
:ord will be used to identify the atom.

% `dataType`
:  Type of information requested:
         Atomic symbol.
 us      Atomic radius.
 r       Color of the atom from the CPK color scheme.
         Average mass of the element.
 name    Name of the element.
         tomic index.
         All atomic data returned in a struct.

### See Also

[sw_mff](sw_mff.html) and [sw_cff](sw_cff.html)

