---
{title: spinw.moment method, link: spinw.moment, summary: calculates the size of the
    reduced moment due to quantum and thermal fluctuations, keywords: sample, sidebar: sw_sidebar,
  permalink: spinw_moment, folder: spinw, mathjax: 'true'}

---

### Syntax

`m = moment(obj,Name,Value)`

### Description

To calculate the reduced moment due to quantum and thermal fluctuations
the magnon poulation is calculated at temperature T over the full
Brillouin zone. To integrate numerically over the Brillouin zone random Q
points are generated and the magnon populations are summed over these
random Q points.
 

### Input Arguments

`obj`
: [spinw](spinw) object.

### Name-Value Pair Arguments

`'nRand'`
: The number of random Q points in the Brillouin-zone,
  default value is 1000.

`'T'`
: Temperature, default value is taken from obj.single_ion.T.

`'tol'`
: Tolerance of the incommensurability of the magnetic
  ordering wavevector. Deviations from integer values of the
  ordering wavevector smaller than the tolerance are
  considered to be commensurate. Default value is 1e-4.

`'omega_tol'`
: Tolerance on the energy difference of degenerate modes when
  diagonalising the quadratic form, default value is 1e-5.

### Output Arguments

'M' is a structure, with the following fields:
moment        Size of the reduced moments, dimensions are [1 nMagExt].
T           	Temperature.
nRand         Number of random Q points.
obj           The copy of the input obj.
Example 1:
tri = sw_model('triAF',1);
M = tri.moment('nRand',1e7);
The example calculates the moment expectation value at zero temperature
on the triangular lattice Heisenberg antiferromagnet. The function gives
the reduced moment (M=0.7385 for S=1 after nRand = 1e7, the reduction is
0.2615). The result can be compared with the following calculations:
[1] A. V Chubukov, S. Sachdev, and T. Senthil, J. Phys. Condens. Matter 6, 8891 (1994).
M = S - 0.261
[2] S. J. Miyake, J. Phys. Soc. Japan 61, 983 (1992).
M = S - 0.2613 + 0.0055/S (1/S is a higher order term neglected here)
Example 2:
sq = sw_model('squareAF',1);
M = sq.moment('nRand',1e7);
The reduced moment of the Heisenberg square lattice antiferromagnet at
zero temperature can be compared to the following published result:
[1] D. A. Huse, Phys. Rev. B 37, 2380 (1988).
M = S - 0.197

### See Also

[spinw](spinw) \| [spinw.spinwave](spinw_spinwave) \| [spinw.genmagstr](spinw_genmagstr) \| [spinw.temperature](spinw_temperature)

{% include links.html %}
