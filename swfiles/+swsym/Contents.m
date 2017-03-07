% The SWSYM library handles the symmetry operations for SpinW objects. The
% functions can be used independently from other parts of SpinW. All
% symmetry operators (symOp) are defined by a matrix with dimensions [3 4
% nOp], where symOp(:,1:3,:) defines the rotation matrices while the
% symOp(:,4,:) the corresponding translations. Also the standard settings
% of the space groups are stored in the symmetry.dat file that can be
% loaded using the SWSYM.generator or SWSYM.operator functions.
%
% Files
%   add       - saves user defined symmetry operators
%   bond      - generates all symmetry equivalent bonds
%   generator - returns symmetry operators of a given space group label or operator string
%   genreduce - reduces the list of symmetry operators to the generators
%   isop      - function determines whether the matrix is a symmetry operator
%   operator  - calculates all symmetry operators or general positions for a space group
%   oporder   - determine the order of the symmetry operator
%   point     - determines point group symmetry at a given position
%   position  - generates symmetry equivalent positions
%   str       - generates a string equivalent of symmetry operators
