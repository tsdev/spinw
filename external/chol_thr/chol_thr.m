% calculates the Cholesky factorisation of a stack of matrices
%
% R     = CHOL_OMP(M);             - Errors if M not positive definite.
% [R,P] = CHOL_OMP(M);             - No error, R is NxNxL;
% L     = CHOL_OMP(M,'lower');
% [L,P] = CHOL_OMP(M,'lower');
%
% The code uses Lapack calls and threads to speed up the calculation. Unlike
% the matlab built-in this mex does not handle sparse matrices.
%
% Input:
%
% M         Matrix with dimensions of [N N L], M=R'*R=L*L'
%
% Options:
%
% tol       Tolerance, if the input is not positive definite, the tol value
%           will be added to each NxN matrix diagonal.
% Colpa     The output is redefined:
%               [K2,invK] = CHOL_OMP(M,'Colpa');
%           where K2 = (R*gComm*R') Hermitian matrix with 
%           gComm = [1...1,-1...1] is the commutator and invK = inv(R).
%
% 
% Note this function does not check for symmetry. Also note that we do not
% truncate the output matrix if the input is not positive definite, unlike
% the built-in. i.e. R or L is always of size (N x N) not (P x P).
% 
% This is a MEX-file for MATLAB.
% 
% Original Author: M. D. Le  [duc.le@stfc.ac.uk]

