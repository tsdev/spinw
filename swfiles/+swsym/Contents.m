% package to handle symmetry operations 
%
% This package deals with symmetry operators of crystallographic space
% groups. It can read the standard space group definitions stored in
% [symmetry.dat], generate all symmmetry elements, determine all symmetry
% equivalent positions, etc. 
%
% All symmetry operators `symOp` are defined by a matrix with dimensions of
% $[3\times 4\times n_{op}]$, where `symOp(1:3,1:3,:)` stores the $[3\times
% 3]$ rotation matrices while the `symOp(1:3,4,:)` holds the corresponding
% translation vectors.
%
% ### Files
%
%   swsym.add      
%   swsym.bond     
%   swsym.generator
%   swsym.genreduce
%   swsym.isop     
%   swsym.operator 
%   swsym.oporder  
%   swsym.point    
%   swsym.position 
%   swsym.str      
