---
{title: spinw.fourier( ), summary: calculates the Fourier transformation of the Hamiltonian,
  keywords: sample, sidebar: sw_sidebar, permalink: spinw_fourier.html, folder: spinw,
  mathjax: 'true'}

---
 
res = FOURIER(obj,hkl,'option1', value1 ...)
 
The function calculates the following sum:
      J(k) = sum_ij J_ij * exp(i*k*d_ij)
The code is optimised for calculating the sum for large number of
k-vectors and alternatively for a large number of d_ij. The single ion
anisotropy is not included in the sum.
 
Input:
 
obj           Input structure, spinw class object.
hkl           Defines the Q points where the Fourier transform is
              calculated, in reciprocal lattice units, size is [3 nHkl].
              Q can be also defined by several linear scan in reciprocal
              space. In this case hkl is cell type, where each element of
              the cell defines a point in Q space. Linear scans are
              assumed between consecutive points. Also the number of Q
              points can be specified as a last element, it is 100 by
              defaults. For example: hkl = {[0 0 0] [1 0 0]  50}, defines
              a scan along (h,0,0) from 0 to 1 and 50 Q points are
              calculated along the scan.
 
              For symbolic calculation at a general reciprocal space
              point use sym class input. For example to calculate the
              spectrum along (h,0,0): hkl = [sym('h') 0 0]. To do
              calculation at a specific point do for example sym([0 1
              0]), to calculate the spectrum at (0,1,0).
 
Options:
 
extend        If true, the Fourier transform will be calculated on the
              magnetic supercell, if false the crystallographic cell will
              be considered. Default is true.
isomode       Defines how Heisenberg/non-Heisenberg Hamiltonians are
              treated. Can have the following values:
                  'off'   Always output the (3x3) form of the
                          Hamiltonian, (default).
                  'auto'  If the Hamiltonian is Heisenberg, only output
                          one of the diagonal values from the (3x3)
                          matrices to reduce memory consumption.
fid           Defines whether to provide text output. Default is defined
              by the swpref.getpref('fid') command. The possible values:
                  0       No text output is generated.
                  1       Text output in the MATLAB Command Window.
                  fid     File ID provided by the fopen() command, the
                          output is written into the opened file stream.
 
Output:
 
res           Structure with the following fields:
  ft          contains the Fourier transfor in a matrix with dimensions
              [3,3,nMagExt,nMagExt,nHKL] or [1,1,nMagExt,nMagExt,nHKL]
              for Heisenberg and non-Heisenberg Hamiltonians respectively
              (if isomode is 'auto'). Here nMagExt is the number of
              magnetic atoms in the magnetic cell and nHKL is the number
              of reciprocal space points.
  hkl         Matrix with the given reciprocal space points stored in a
              matrix with dimensions [3,nHKL].
  isiso       True is the output is in Heisenberg mode, when the ft
              matrix has dimensions of [1,1,nMagExt,nMagExt,nHKL],
              otherwise is false.
 
See also SPINW.OPTMAGK.
 

