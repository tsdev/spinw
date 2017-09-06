---
{title: addcoupling( ), keywords: sample, summary: assigns a predefined matrix as exchange coupling on selected bonds,
  sidebar: sw_sidebar, permalink: '@spinw_addcoupling.html', folder: '@spinw', mathjax: 'true'}

---
  assigns a predefined matrix as exchange coupling on selected bonds
 
  ADDCOUPLING(obj, 'option1', value1, ...)
 
  Input:
 
  obj           spinw class object.
 
  Options:
 
  mat           Label or index of the matrix that will be assigned to
                selected bonds.
  bond          Selects the interacting atom pairs through the
                obj.coupling.idx number. The coupling.idx numbers are in
                increasing order according to the distances between
                magnetic atoms, for example all shortest interatom
                distances have idx=1, second shortest idx=2 and so on.
                bondIdx can be a vector to assign the matrix to multiple
                inequivalent bonds.
  atom          Contains labels of atoms or index of atoms in the
                obj.unit_cell list of atoms. If a single string label or
                number is given, only bonds between the selected atoms will
                be assigned. If a cell with 2 strings are given, only bonds
                between the two selected atoms will be assigned. Only works
                if space group is P0. Default is [].
  subIdx        If the above options are not enough to select the desired
                bonds, using subIdx bonds can be selected one-by-one from
                the list of bonds that fulfill the above options. Only
                works if the space group is P0.
  type          Type of the coupling. Possible values:
                    0       Quadratic exchange, default.
                    1       Biquadratic exchange.
                Can be also one of the following string: 'quadratic',
                'biquadratic'.
  sym           If true, symmetry operators will be applied on the exchange
                matrices to generate the coupling on symmetry equivalent
                bonds, if false all symmetry equivalent bonds will have the
                same exhcange matrix. It has to be false if subIdx is
                given.
 
  Output:
 
  The function adds extra entries in the 'coupling.matrix' field of the obj
  spinw object.
 
  Example:
 
  ...
  cryst.addmatrix('label','J1','value',0.123)
  cryst.gencoupling
  cryst.addcoupling('mat','J1','bond',1)
 
  This will add the 'J1' diagonal matrix to all second shortes bonds
  between magnetic atoms.
 
  See also SPINW, SPINW.GENCOUPLING, SPINW.ADDMATRIX.
 
