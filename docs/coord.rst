Coordinate systems
==================

SpinW uses three different coordinate systems. It is important to consider every time, that a given vector/tensor data is in which coordinate system. In the following documentation the xyz, lu and rlu shorthand notation will be used for the three coordinate systems.

Lattice units (lu) coordinate system
-------------------------------------

This is the lattice coordinate system. The three axis are the **a**, **b** and **c** crystal axes. For example the vector :code:`[0,1,0]` denotes a position translated from the origin by **b** vector. The following sw class properties are stored in lattice units:

* atomic positions (:attr:`sw.unit_cell`.r)
* translation vectors for bonds (:attr:`sw.coupling`.dl)

Also several function takes input or output in lattice units:

* atomic positions of the output of :func:`matom` and :func:`atom` methods of sw
* magnetic moments can be given in lattice units for the :func:`genmagstr` function (using the 'unitS' option with 'lu' value)
* calculated bond vectors from :func:`couplingtable`

xyz coordinate system
---------------------

Most of the sw class properties are stored in the xyz coordinate system. The xyz coordinate system is right-handed Cartesian and fixed to the
crystal lattice:

* *x*: parallel to *a*-axis,
* *y*: perpendicular to *x* and in the *ab*-plane,
* *z*: perpendicular to *xy*-plane

The following properties are in the xyz coordinate system:

* twin rotation matrices (:attr:`sw.twin`.rotc)
* stored 3x3 matrices (:attr:`sw.matrix`.mat)
* magnetic field (:attr:`sw.single_ion`.field)
* magnetic moment components (:attr:`sw.mag_str`.S)
* normal vector of the magnetic structure (:attr:`sw.mag_str`.n)

Also output of several functions are in xyz coordinate system:

* spin-spin correlation function matrix elements calculated by :func:`spinwave` method (spec.Sab matrices)
* interaction matrices calculated by :func:`couplingtable`

Vectors in reciprocal space in |AA|:math:`^{-1}` units are also in a right handed Cartesian coordinate system:

* momentum transfer values in |AA|:math:`^{-1}` reciprocal lattice units of the spec.hklA output of the :func:`spinwave` function

Reciprocal lattice (rlu) coordinate system
------------------------------------------

The reciprocal lattice coordinate system is the dual vector space of the
lattice coordinate system. The three axis are the reciprocal lattice vectors denoted by **a**:math:`^*`, **b**:math:`^*` and **c**:math:`^*`. The following sw class properties are stored in rlu units:

* magnetic ordering wave vector (:attr:`sw.mag_str`.k)

Also several function takes input in rlu units:

* the :func:`sw_neutron` function takes the option 'uv', that defines the scattering plane by two vectors in rlu
* the first input of the :func:`spinwave` function is a list of Q points in rlu units

Transformation between coordinate systems
-----------------------------------------

To transform vectors and tensors between the above coordinate systems, the output of
:func:`basisvector` function can be used::

	tri = sw;
	tri.genlattice('lat_const',[3 3 5],'angled',[90 90 120])
	BV = tri.basisvector

The :code:`BV` 3x3 matrix contains the coordinates of the lattice vectors as column vectors: :code:`[a b c]`. To convert vectors from the abc coordinate system to the xyz, we can do the following::

	r_abc = [1/2 1/2; 0];
	r_xyz = BV * r_abc

.. |AA| unicode:: U+212B
