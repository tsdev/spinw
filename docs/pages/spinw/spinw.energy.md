---
{title: spinw.energy( ), summary: calculates the ground state energy per spin, keywords: sample,
  sidebar: sw_sidebar, permalink: spinw_energy.html, folder: spinw, mathjax: 'true'}

---
calculates the ground state energy per spin
 
E = ENERGY(obj, Option1, Value1 ...)
 
The extended magnetic unit cell, stored in obj, is used for the
calculation. For non-zero k vector, the interaction energies between
neighbouring extended unit cells depend on the direction of the moments
in the two extended unit cells. The angles in further extended unit cells
are calculated based on the k vector (the k vector is in the units of the
crystallographic unit cell) and the n vector (normal to the spin rotation
plane). The moment directions in further extended unit cells are
calculated by rotating the spins of the extended unit cell in the origin
by k*R degree around the n vector, where R is the translation vector of
the origin of the farther extended unit cell. If the extended unit cell
is equivalent to the crystallographic unit cell, this is equivalent to
the standard definition of the k vector.
 
Input:
 
obj       spinw class object.
 
Options:
 
epsilon   The smallest value of incommensurability that is tolerated 
          without warning. Default is 1e-5.
 
Output:
 
E         Energy per moment (anisotropy + exchange + Zeeman energy).
 
 
WARNING WARNING WARNING WARNING WARNING WARNING WARNING WARNING WARNING
 
The calculated energy can be wrong for incommensurate structures. For
example a structure where the spins are rotating in XY plane with an
incommensurate wavevector of (1/3,0,0). The function only calculates the
anisotropy energy in the first unit cell, that is for single spin
Eaniso = Axx*Sxx^2+Ayy*Syy^2. While the anisotropy energy in reality is
independent of the spin orientation in the XY plane Eaniso=3S*(Axx+Ayy)/2.
Thus for incommensurate structures one has to be carefull! In the
triangular case one has to extend the unit cell to nExt = [3 3 1] (in the
hexagonal setting), in this case the energy will be correct.
 
WARNING WARNING WARNING WARNING WARNING WARNING WARNING WARNING WARNING
 
Example:
 
...
cryst.optmagstr('nRun',10)
E = cryst.energy
 
After optimising the magnetic structure (by minimizing the ground state 
energy), the energy per spin is calculated. This can be compared to
different ground state structures to decide which is the right classical
ground state of the magnetic model in cryst.
 
See also SPINW, SPINW.ANNEAL, SPINW.NEWCELL.
 
