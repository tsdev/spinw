%% Coordinate systems
% SpinW uses three different coordinate systems.

%% abc coordinate system
% This is the lattice coordinate system, every vector, whose components are
% given in lattice units are in this coordinate system. The three axis are
% the *a*, *b* and *c* crystal axes. The axis length can be normalized to
% one or to the lattice parameters. The following sw class properties are
% stored in lattice units:
%
% * atomic positions (*sw.unit_cell.r*)
% * translation vectors for bonds (*sw.coupling.dl*)
%
% Also several function takes input or aoutput in lattice units:
%
% * atomic positions of the output of sw.matom and sw.atom methods
%   (sw.matom.r, sw.atom.r)
% * magnetic moments can be given in lattice units for the sw.genmagstr
%   method (using the 'unitS' option with 'lu' value)
% * calculated bond vector by sw.couplingtable
%

%% xyz coordinate system
% Most of the sw class properties are stored in the xyz coordinate system.
% The xyz coordinate system is right-handed Cartesian and fixed to the
% crystal lattice:
%
% * *x*: parallel to *a*-axis,
% * *y*: perpendicular to *x* and in the *ab*-plane,
% * *z*: perpendicular to *xy*-plane
%
% The following properties are in xyz coordinate system:
%
% * twin rotation matrices (*sw.twin.rotc*)
% * stored 3x3 matrices (*sw.matrix.mat*)
% * magnetic field (*sw.single_ion.field*)
% * magnetic moment components (*sw.mag_str.S*)
% * normal vector of the magnetic structure (*sw.mag_str.n*)
%
% Also output of different functions are in xyz coordinate system:
%
% * spin-spin correlation function matrix elements calculated by
%   sw.spinwave method (spec.Sab matrices)
% * interaction matrices calculated by sw.couplingtable
%
% Also any vector in reciprocal space that is in Angstrom^-1 units:
%
% * momentum transfer values in Angstrom^-1 units of the calculated
%   spectrum using sw.spinwave function (spec.hklA)

%% Reciprocal lattice coordinate system
% The reciprocal lattice coordinate system is the dual vector space of the
% lattice coordinate system. The three axis are the reciprocal lattice
% vectors denoted by *a**, *b** and *c**. The following sw properties are
% stored in r.l.u.:
%
% * magnetic ordering wave vector (sw.mag_str.k)
%
% Also several function takes input in r.l.u.:
%
% * the sw_neutron function takes the option 'uv', that defines the
%   scattering plane by two vectors in r.l.u.
% * the first input of the sw.spinwave function is a list of Q points in
%   r.l. units
% * the sw.genmagstr function can take magnetic moment components in r.l.
%   units if the 'unitS' option is set to 'lu'
%

%% Transformation between coordinate systems
% To transform between the above coordinate systems, the output of
% sw.basisvector can be used:

tri = sw;
tri.genlattice('lat_const',[3 4 5],'angled',[90 90 120])
BV = tri.basisvector

%%
% The BV 3x3 matrix contains the lattice vectors as column vectors: [a b
% c]. To convert from abc coordinate system to the xyz, we need 3x1 column
% vectors:

r_abc = [1/2 1/2; 0];
r_xyz = BV * r_abc

%%
% 








