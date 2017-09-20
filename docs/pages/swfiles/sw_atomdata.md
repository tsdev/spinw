---
{title: sw_atomdata( ), link: sw_atomdata, summary: returns information on elements
    stored in the atom.dat file, keywords: sample, sidebar: sw_sidebar, permalink: sw_atomdata.html,
  folder: swfiles, mathjax: 'true'}

---

### Syntax

`[data atomlabel] = sw_atomdata(atomsymb, {datatype})`

### Description



### Examples

sw_atomdata('H','radius') = 0.37
If the atom label does not exists, the function returns radius = 1,
color = [255 167 0].
optional second output is 'atomLabel' that contains the name of the atom
clean.

### Input Arguments

`atomSymb`
: String of the name of the atom, for example 'He' or the atomic
  number Z. If the string contains whitespace character, the
  second word will be used to identify the atom.

`dataType`
: Type of information requested:
      name        Atomic symbol.
      radius      Atomic radius.
      color       Color of the atom from the CPK color scheme.
      mass        Average mass of the element.
      longname    Name of the element.
      Z           tomic index.
      all         All atomic data returned in a struct.

### See Also

[sw_mff](sw_mff.html) \| [sw_cff](sw_cff.html)

