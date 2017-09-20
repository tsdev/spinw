---
{title: sw_readspec, link: sw_readspec, summary: read spin wave dispersion data from
    file, keywords: sample, sidebar: sw_sidebar, permalink: sw_readspec.html, folder: swfiles,
  mathjax: 'true'}

---

### Syntax

`data = sw_readspec(datapath)`

### Description

It reads experimental spin wave dispersion data from text file, whose
location is defined in path.
 
Format of the input data file:
Every line of the data file contains information about an energy scan at
constant Q. Data consists of floating point numbers separated by space
(first line can be a header line):
 
  QH QK QL minE maxE I1 E1 s1 I2 E2 s2 ...
      where:
 
QH        H index of the Q point,
QK        K index of the Q point,
QL        L index of the Q point,
minE      lower boundary of the E scan,
maxE      upper boundary of the E scan,
In        intensity of the n-th spin wave mode,
En        center of the n-th spin wave mode, has to be in increasing
          order,
sn        standard deviation of the corresponding energy
 
The number of modes in a single line of the data file is unlimited,
however in every line the number of modes have to be the same. Scans with
less modes should contain zero intensities. The code automatically omits
modes that have either zero intensity or zero energy.
 
Before any data line a special line can be inserted that gives the
measured correlation in square brackets, for axample: [Mxx+Myy]. For the
formatting of this string, see <a href="matlab:doc sw_parstr">sw_parstr</a>.
If the measured type of correlation is undefined, unpolarised neutron
scattering intensity is assumed ([Sperp]). When cross sections measured
in the Blume-Maleev coordinate system (see <a href="matlab:doc sw_egrid">sw_egrid</a>), the normal to the
scattering plane has to be also defined. This can be given in a second
pair of square brackes in the xyz coordinate system, for example: [Myy]
[1 0 0]. If n is undefined, the default value is [0 0 1].
 
 
Example input data file (polarised scans in the (0KL) plane):
 
QH    QK        QL      ENlim1  ENlim2  I1  EN1       s1    I2  EN2       s2
[Mxx] [1 0 0]
0     1        2.9992   0       15      1    3.7128   1.0   1   8.6778    1.0
0     1        2.8993   0       15      1    7.0000   1.0   1   11.1249   1.0
0     1        2.7993   0       20      1   13.8576   1.0   0   0.0       0.0
0     1        2.6994   0       20      1   17.3861   1.0   0   0.0       0.0
[Myy] [1 0 0]
0     1.0000   2.0000   0       25      1   20.2183   1.0   0   0.0       0.0
0     1.1000   2.0000   15      30      1   22.7032   1.0   0   0.0       0.0
0     1.2000   2.0000   20      35      1   25.1516   1.0   0   0.0       0.0
 

### See Also

[sw_egrid](sw_egrid.html) \| [spinw.fitspec](spinw_fitspec.html)

