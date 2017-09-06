---
{title: addatom( ), keywords: sample, summary: adds new atom to an spinw object, sidebar: sw_sidebar,
  permalink: '@spinw_addatom.html', folder: '@spinw', mathjax: 'true'}

---
  adds new atom to an spinw object
 
  ADDATOM(obj,'Option1', Value1, ...)
 
  Input:
 
  obj       spinw object
 
  Options:
 
  r         Atomic positions, dimensions are [3 nAtom]. No default value!
  S         Spin of the atoms, dimensions are [1 nAtom], for non-magnetic
            atoms set S to zero. Default spin is generated from the given
            label of the atom. For example if 'label' is 'MCr3+' or 'Cr3+'
            then the high spin of S=3/2 is automatically generated. The
            high spin values for every ion is stored in the last column of
            the ion.dat file. If the atom type is unknown S=0 is assumed.
  label     Names of the atoms for plotting and form factor
            calculations (see ion.dat), it is a cell, optional.
            Example:
            {'atom1' 'atom2' 'atom3'}
            Default value is 'atomI', where I is the atom index.
  color     Colors of the atoms for plotting, dimensions are [3 nAtom],
            where each column describes an RGB color. Each value is between
            0 and 255, optional. Default value is [255;165;0] for each
            atom.
            Alternatively a name of the color can be given as a string, for
            example 'White', for multiple atoms package it into a cell. For
            the list of colors, see swplot.color().
  ox        Oxidation number given as a double or it will be determined
            automatically from label. Default is 0.
  occ       Occupancy, given as double. Default is 1.
  formfact  Neutron scattering form factor, given as 9 numbers, for details
            see the help of sw_mff(). Also string labels can be used.
  formfactn  Neutron scattering form factor, given as 9 numbers, for details
            see the help of sw_mff(). Also string labels can be used.
  formfactx X-ray scattering form factor, given as 9 numbers, for details
            see the help of sw_cff().
  Z         Atomic number, given as integer or determined from label
            automatically. Default is 113 (Unobtanium).
  A         Atomic mass, given as integer. Default is -1 for the natural
            mixture of isotopes.
  bn        Neutron scattering length, given as double.
  bx        X-ray scattering length.
  biso      Isotropic displacement factors in units of Angstrom^2.
            Definition is the same as in FullProf, defining the
            Debye-Waller factor as:
                Wd = 1/8*biso/d^2
            including in the structure factor as exp(-2Wd).
  update    If true, existing atom with the same label and position as a
            new one will be updated. Default is true.
    
  Output:
 
  The function creates extra elements in the 'unit_cell' field of the obj
  spinw object.
 
  Example:
 
  crystal.ADDATOM('r',[0 1/2; 0 0; 0 0],'S',[1 0],'color',{'red' 'blue'})
 
  Adds a magnetic atom (S=1) at position (0,0,0) and a non-magnetic one at
  (1/2 0 0) with red and blue color respectively.
 
  See also SPINW.GENLATTICE, SPINW.ADDMATRIX, SWPLOT.COLOR, SW_MFF, SW_CFF.
 
