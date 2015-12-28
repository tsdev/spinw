% Diagonalises a stack of matrices
%
% [V, D] = EIG_THR(H,'orth','sort')
%
% The code uses Lapack calls and pthreads/windows threads to speed up the 
% calculation.
%
% Input:
%
% H     Matrix with dimensions [N N L].
%
% Options:
%
% orth  Orthogonalize the eigenvectors in case of degenerate eigenvalues
%       for non-Hermitian matrices (for Hermitian matrices the eigenvectors
%       are always orthogonal).
% sort  Sort the eigenvalues in increasing order (and sort the
% corresponding eigenvectors as well).
%
% Output:
%
% V     Matrix with dimesions [N N L], contains the eigenvectors for each
%       NxN matrix in columns.
% D     Matrix with dimensions [N L], each column contain the corresponding
%       eigenvalues.
%  
% This is a MEX-file for MATLAB.
%  
% In addition the the mexFunction, the following subsidiary functions are
% defined:
%  
%   void fliplr()    - Flips matrix by columns
%   void quicksort() - Quicksort algorithm, returns indices
%   void sort()      - Sort and permutes eigenvalue/vectors
%   void orth()      - Orthogonalises degenerate eigenvecs
%   
% Original Author: M. D. Le  [duc.le@stfc.ac.uk]

