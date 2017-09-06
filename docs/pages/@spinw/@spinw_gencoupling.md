---
{title: gencoupling( ), keywords: sample, summary: generates the COUPLING property of spinw object,
  sidebar: sw_sidebar, permalink: '@spinw_gencoupling.html', folder: '@spinw', mathjax: 'true'}

---
  generates the COUPLING property of spinw object
 
  GENCOUPLING(obj, 'option1', value, ...)
 
  The function calculates equivalent bonds between magnetic atoms. These
  are determined either based on crystal symmetry or bond length (with
  tolDist tolerance). If the space group index of 0 is defined
  (obj.lattice.sym=0), the equivalent bonds will be based on bond length.
  For space group index larger than 0, the symmetry equivalent bonds will
  be determined. This can ve overwritten by the forceNoSym parameter to
  consider bond length.
 
  IMPORTANT!
    This function has to be used after the crystal structure is defined.
    The SPINW.ADDCOUPLING, SPINW.COUPLINGTABLE functions will only work
    afterwards.
 
  Input:
 
  obj           spinw class object.
 
  Options:
 
  forceNoSym    If true, equivalent bonds are generated based on
                bond length with .tolDist tolerance. If false symmetry
                operators will be used if they are given
                (obj.lattice.sym>0).
  maxDistance   Maximum bond length that will be stored in the
                obj.coupling property in units of Angstrom. Default is 8.
  maxSym        Maximum bond length until the symmetry equivalent bonds are
                generated. It is usefull if long bonds have to be generated
                for the dipolar interaction, but the symmetry analysis of
                them is not necessary. Default value is equal to
                maxDistance.
  tolDist       Tolerance of distance, within two bonds are regarded
                equivalent, default is 1e-3 Angstrom. Only used, when no
                space group is defined.
  dMin          Minimum bond length, below which an error is triggered.
                Default value is 0.5 Angstrom.
 
  Output:
 
  The obj.coupling field will be filled with values, depending on the
  crystal geometry.
 
  Example:
 
  cryst = spinw;
  cryst.genlattice('lat_const',[3 3 5],'angled',[90 90 120])
  cryst.addatom('r',[0 0 0])
  cryst.gencoupling
  cryst.couplingtable(1:3)
 
  A triangular lattice is created in cryst and after using gencoupling()
  the couplingtable() function lists the 1st, 2nd and 3rd neighbor bonds.
 
  See also SPINW, SPINW.SYMMETRY, SPINW.NOSYM.
 
