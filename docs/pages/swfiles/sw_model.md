---
{title: sw_model, link: sw_model, summary: creates different predefined spin models,
  keywords: sample, sidebar: sw_sidebar, permalink: sw_model.html, folder: swfiles,
  mathjax: 'true'}

---

### Syntax

`obj = sw_model(model, param, {fid})`

### Description



### Input Arguments

`model`
: String, name of the model, one of the following:
      'triAF'     Triangular lattice Heisenberg antiferromagnet
                  in the ab plane (a=b=3 Å), with γ =
                  120 deg angle and optimised magnetic structure.
                  Arbitrary number of Heisenberg interaction can
                  be defined, param(1) gives the value of 1st
                  neighbor interaction, param(2) the second etc.
      'squareAF'  Square lattice antiferromagnet.
      'chain'     Chain with further neighbor interactions.

`param`
: Input parameters of the model, depending on which is selected.

`fid`
: Where to print the text output. Default is 1 to print to the
  Command Window. Optional.

### Output Arguments

obj       spinw class object with the selected model.

### See Also

[spinw](spinw.html)

