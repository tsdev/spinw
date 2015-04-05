%% SW Class Properties
% The sw object properties store all the information necessary for the spin
% wave calculation. It has 8 fields with several subfields, see below.

%% lattice
% The *lattice* field stores the crystallographic unit cell parameters.
% Subfields are:
%
% * *lat_const*: Lattice constants in a vector with dimensions of [1 3] in
%               Angstrom units.
% * *angle*:    Alpha, beta, gamma angles in a vector with dimensions of
%               [1 3] in radian units.
% * *sym*:      Crystal space group, integer number that denotes the line
%               number in the symmetry.dat file
%
% See also SW.GENLATTICE, SW.ABC, SW.BASISVECTOR, SW.NOSYM.

%% unit_cell
% The *unit_cell* field stores the information on the atoms in the
% crystallographic unit cell. Subfields are:
%
% * *r*:        Atomic positions, matrix with dimensions of [3 nAtom], in
%               lattice units.
% * *S*:        Spin quantum number of the atoms, vector with dimensions of
%               [1 nAtom], non-magnetic atoms have S=0, can be also the J
%               quantum number.
% * *label*:    Label of the atoms, strings in a cell with dimensions of [1
%               nAtom].
% * *color*:    Color of the atoms for plotting in a matrix with dimensions
%               of [3 nAtom], every column is an RGB code (numbers between
%               0-255).
%
% See also SW.ADDATOM, SW.ATOM, SW.MATOM, SW.NEWCELL, SW.PLOT.

%% twin
% The *twin* field stores information on crystallographic twins. The
% subfields are:
%
% * *rotc*:     Rotation matrices in the xyz coordinate system for
%               every twin, stored in a matrix with dimensions of [3 3
%               nTwin].
% * *vol*:      Volume ratio of the different twins, stored in arow vector
%               with nTwin elements.
%
% See also SW.ADDTWIN, SW.TWINQ, SW.UNIT_CELL.

%% matrix
% The *matrix* field stores 3x3 matrices that can be assigned to the
% magnetic Hailtonian. The subfields are:
%
% * *mat*:      It stores the actual values of 3x3 matricesstacked along the
%               3rd dimension.
% * *color*:    RGB color codes assigned for every matrix, stored in a
%               matrix with dimensions of [3 nMatrix], each column is an
%               [R;G;B] codewith numbers between 0 and 255.
% * *label*:    Label for every matrix, stored as strings in a
%               cell with dimensions of [1 nCell].
%
% See also SW.ADDMATRIX, SW.NTWIN.

%% single_ion
% The *single_ion* field stores single ion terms of the Hamiltonian.
% The subfields are:
%
% * *aniso*:    Row vector contains nMagAtom integers, each integer
%               assignes one of the matrices from the sw.matrix field
%               to a magnetic atom in the sw.matom list as a single
%               ion anisotropy (zeros for no anisotropy).
% * *g*:        Row vector of nMagAtom integers, each integer
%               assignes one of the matrices from the sw.matrix field
%               to a magnetic atom in the sw.matom list as a
%               g-tensor
% * *field*:    External magnetic field stored in a row vector with 3
%               components in the xyz coordinate system, default unit is
%               Tesla.
% * *T*:        Temperature, scalar, default unit is Kelvin.
%
% See also SW.ADDANISO, SW.ADDG, SW.GETMATRIX, SW.SETMATRIX, SW.INTMATRIX.

%% coupling
% The *coupling* field stores the list of bonds, where the bond is
% between atom1 and atom2. The subfields are:
%
% * *dl*:       The distance between the unit cells of two interacting
%               spins, stored in a matrix with dimensions of [3 nCoupling].
% * *atom1*:    First magnetic atom, pointing to the list of
%               magnetic atoms in sw.matom list, stored in a row vector
%               with nCoupling elements.
% * *atom2*:    Second magnetic atom, with same dimensions.
% * *mat_idx*:  Integers selecting  exchange matrices from sw.matrix field
%               for every bond in stored in a matrix with dimensions of [3
%               nCoupling], maximum three matrices per coupling (zeros for
%               no coupling).
% * *idx*:      Increasing indices for the symmetry equivalent
%               couplings, starting with 1,2,3 (equivalent to first,
%               second, third... neighbor).
%
% See also SW.GENCOUPLING, SW.ADDCOUPLING, SW.FIELD.

%% mag_str
% The *mag_str* field stores the magnetic structure. The subfields are:
%
% * *S*:        It stores the moment direction for every magnetic atom of the
%               magnetic supercell (can be equivalent to the
%               crystallographic unit cell) in a matrix with dimensions of [3 nMagExt],
%               where nMagExt = nMagAtom*prod(mag_str.N_ext). The moment
%               Components are in the xyz coordinate system.
% * *k*:        Magnetic ordering wave vector stored in a row vector with 3
%               components in reciprocal lattice units.
% * *n*:        Normal vector to the rotation of the moments in
%               case of non-zero ordering wave vector (helical magnetic
%               structures), row vector with three elements. Components are
%               in the xyz coordinate system.
% * *N_ext*:    Size of the magnetic supercell, default is [1 1 1]
%               if the magnetic cell is identical to the
%               crystallographic cell, the 1x3 vector extends the
%               cell along the a, b and c axis
%
% See also SW.GENMAGSTR, SW.OPTMAGSTR, SW.ANNEAL, SW.MOMENT, SW.NMAGEXT, SW.STRUCTFACT.

%% unit
% The *unit* field stores the conversion factors between energy, magnetic
% field and temeprature units in the Hamiltonian. Defaults units are meV,
% Tesla and Kelvin for energy, magnetic field and temperature. The
% subfields are:
%
% * *kB*:       Boltzmann constant, default is 0.0862 (meV/K).
% * *muB*:      Bohr magneton, default is 0.0579 (meV/T).
%