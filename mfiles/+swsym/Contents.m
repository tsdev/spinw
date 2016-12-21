% The SWSYM library handles the symmetry operations for SpinW objects. The
% functions can be used independently from other parts of SpinW. All
% symmetry operators (symOp) are defined by a matrix with dimensions [3 4
% nOp], where symOp(:,1:3,:) defines the rotation matrices while the
% symOp(:,4,:) the corresponding translations. Also the standard settings
% of the space groups are stored in the symmetry.dat file that can be
% loaded using the SWSYM.generator or SWSYM.operator functions.
%
% List of functions:
%
% <a href="matlab: doc swsym.add">add</a>           Saves user defined symmetry operators.
% <a href="matlab: doc swsym.bond">bond</a>          Generates all symmetry equivalent bonds.
% <a href="matlab: doc swsym.generator">generator</a>     Returns symmetry operators of a given space group.
% <a href="matlab: doc swsym.genreduce">genreduce</a>     Reduces the list of symmetry operators to the generators.
% <a href="matlab: doc swsym.isop">isop</a>          Function determines whether the matrix is a symmetry operator.
% <a href="matlab: doc swsym.operator">operator</a>      Calculates all symmetry operators or general positions for a space group.
% <a href="matlab: doc swsym.oporder">oporder</a>       Determines the order of the symmetry operator.
% <a href="matlab: doc swsym.point">point</a>         Determines point group symmetry at a given position.
% <a href="matlab: doc swsym.position">position</a>      Generates symmetry equivalent positions.
% <a href="matlab: doc swsym.str">str</a>           Generates a string equivalent of symmetry operators.
%
