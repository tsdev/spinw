%% Coordinate systems
% SpinW uses four coordinate systems.

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
% * spin-spin correlation function calculated by sw.spinwave method
%   (spec.Sab matrices)
% * interaction matrices calculated by sw.couplingtable
%

%% reciprocal lattice coordinate system
% The reciprocal lattice coordinate system is the dual vector space of the
% lattice coordinate system. The three axis are the reciprocal lattice
% vectors denoted by *a**, *b** and *c** LATEXf(x,y)=x^2+y^2PATEX.






