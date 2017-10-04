---
{title: spinw.disp method, link: spinw.disp, summary: prints information, keywords: sample,
  sidebar: sw_sidebar, permalink: spinw_disp, folder: spinw, mathjax: 'true'}

---
  
### Syntax
  
`{swdescr} = disp(obj)`
  
### Description
  
`{swdescr} = disp(obj)` generates text summary of a [spinw](spinw) object.
Calling it with output argument, it will generate a text version of the
internal data structure giving also the dimensions of the different
matrices.
  
### Examples
  
Here the internal data structure is generated:
 
```matlab
crystal = spinw
swFields = disp(crystal)
```
*Output*
```
swFields =
    'spinw object (symbolic: off, symmetry: off, textoutput: "stdout")
     lattice
              angle: [1x3 double]
          lat_const: [1x3 double]
                sym: [3x4xnSymOp double]  nSymOp=0
             origin: [1x3 double]
              label: [1xnStr char]
     unit_cell
                  r: [3xnAtom double]  nAtom=0
                  S: [1xnAtom double]
              label: [1xnAtom char]
              color: [3xnAtom integer]
                 ox: [1xnAtom double]
                occ: [1xnAtom double]
                  b: [2xnAtom double]
                 ff: [2x11xnAtom double]
                  A: [1xnAtom integer]
                  Z: [1xnAtom integer]
               biso: [1xnAtom double]
     twin
                vol: [1xnTwin double]  nTwin=1
               rotc: [3x3xnTwin double]
     matrix
                mat: [3x3xnMat double]  nMat=0
              color: [3xnMat integer]
              label: [1xnMat char]
     single_ion
              aniso: [1xnMagAtom integer]  nMagAtom=0
                  g: [1xnMagAtom integer]
              field: [1x3 double]
                  T: [1x1 double]
     coupling
                 dl: [3xnBond integer]  nBond=0
              atom1: [1xnBond integer]
              atom2: [1xnBond integer]
            mat_idx: [3xnBond integer]
                idx: [1xnBond integer]
               type: [3xnBond integer]
                sym: [3xnBond integer]
               rdip: [1x1 double]
               nsym: [1x1 integer]
     mag_str
                  F: [3xnMagExtxnK double]  nK=0
                  k: [3xnK double]
               nExt: [1x3 integer]
     unit
                 kB: [1x1 double]
                muB: [1x1 double]
                mu0: [1x1 double]
              label: [1x4 char]
           nformula: [1x1 integer]
               qmat: [3x3 double]
     '
```
 
  
### Input Arguments
  
`obj`
: [spinw](spinw) object.
  
### Output Arguments
  
`swdescr`
: If output variable is given, the description of the `obj` object
  will be output into the `swdescr` variable, instead of being
  written onto the Command Window/file. Optional.
  
### See Also
  
[spinw](spinw)
 

{% include links.html %}
