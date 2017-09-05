---
title: intmatrix( )
keywords: sample
summary: "creates the interactions matrices (connectors and values)"
sidebar: product1_sidebar
permalink: intmatrix.html
folder: @spinw
mathjax: true
---
  creates the interactions matrices (connectors and values)
 
  [SS, SI, RR] = INTMATRIX(obj, 'Option1', Value1, ...)
 
  Input:
 
  obj           Input spinw class object.
 
  Options:
 
  fitmode       Can be used to speed up calculation, modes:
                    true    No speedup (default).
                    false   For the interactions stored in SS, only the
                            'all' field is calculated.
  plotmode      If true, additional rows are added to SS.all, to identify
                the couplings for plotting. Default is false.
  sortDM        If true each coupling is sorted for consistent plotting of
                the DM interaction. Sorting is based on the dR distance
                vector, pointing from atom1 to atom2. Its components should
                fulfill the following rules in hierarchical order:
                    1. dR(x) > 0
                    2. dR(y) > 0
                    3. dR(z) > 0.
                Default is false.
  zeroC         Whether to output bonds with assigned matrices that are
                zero. Default is false.
  extend        If true, all bonds in the magnetic supercell will be
                generated, if false, only the bonds in the crystallographic
                unit cell is calculated. Default is true.
  conjugate     Introduce the conjugate of the couplings (atom1 and atom2
                exchanged). Default is false.
 
  Output:
 
  SS            Structure with  fields {iso,ani,dm,gen,bq,dip}. It
                describes the interactions between spins. Every field is a
                matrix, where every column is a coupling between two spins.
                The first 3 rows contain the unit cell translation vector
                between the interacting spins, the 4th and 5th row contains
                the indices of the two interacting spins in the 'spin'
                variable. The following rows contains the strength of the
                interaction. For isotropic exchange it is a single number,
                for DM interaction [DMx; DMy; DMz], for anisotropic
                interaction [Jxx; Jyy; Jzz] and for general interaction
                [Jxx; Jxy; Jxz; Jyx; Jyy; Jyz; Jzx; Jzy; Jzz] and for
                biquadratic exchange it is also a single number.
                For example:
                 SS.iso = [dLatX; dLatY; dLatZ; spinIdx1; spinIdx2; Jval].
                For plotmode true, two additional rows are added to SS.all,
                that contains the idx indices of the obj.matrix(:,:,idx)
                corresponding matrix for each coupling and the .idx values
                of the couplings. The dip field contains the dipolar
                interactions that are not added to the SS.all field.
 
  SI            Single ion energy, due to anisotropy and magnetic field.
  SI.aniso      Matrix with dimensions of [3 3 nMagAtom] sized matrix,
                where the energy of the i-th spin is
                E_aniso = spin(:)*A(:,:,i)*spin(:)'.
  SI.g          g-tensor, with dimensions of [3 3 nMagAtom]. It determines
                the energy of the magnetic moment in external field:
                E_field = B(:)*g(:,:,i)*spin(:)'.
  SI.field      External magnetic field [Bx By Bz].
 
  RR            Positions of the atoms in lattice units, dimensions are
                [3 nMAgExt].
 
  See also SPINW.COUPLINGTABLE, SPINW.SYMOP.
 
