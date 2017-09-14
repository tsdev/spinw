% package to handle symmetry operations 
%
% The functions can be used independently from other parts of spinw. All
% symmetry operators `symOp` are defined by a matrix with dimensions of
% $[3\times 4\times n_{op}]$, where `symOp(:,1:3,:)` defines the rotation
% matrices while the `symOp(:,4,:)` the corresponding translations. Also
% the standard settings of the space groups are stored in the
% `symmetry.dat` file that can be loaded using the [swsym.generator] or
% [swsym.operator] functions.
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
