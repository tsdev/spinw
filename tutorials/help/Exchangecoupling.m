%% Exchange coupling
% Exchange couplings define magnetic interactions on selected bonds.
%
%% Definition
% Exchange coupling values similarly to single ion anisotropy matrix and
% g-tensor are represented as a matrix with dimensions of 3x3. Thes value
% of matrices are defined in sw.matrix.mat. These matrices can be assigned
% to arbitrary pair of interacting spins, these assignments are stored in
% sw.coupling.
%
%% Symmetry transformations
% If the exchange interactions are assigned to bonds using symmetry
% operators (the sw.symmetry property is true) the exchange matrices are
% transformed as tensors between symmetry equivalent bonds. The list of
% symmetry equivalent bonds can be acquired using the sw.couplingtable
% command. For example the shortest equivalent bonds (first neigbors) can
% be listed using sw.couplingtable(1).table command:

sq = sw;
sq.genlattice('lat_const',[4 4 5],'sym','P 4')
sq.addatom('r',[0 0 0],'S',1)
sq.gencoupling
sq.couplingtable(1).table

%%
% The above example defines a square lattice of spin-1 magnetic atoms with
% fourfold symmetry in the *ab*-plane. The output of sw.couplingtable gives
% a matrix where each column defines individual couplings. The first three
% rows define the translation between the two interacting spins in lattice
% units. So the first bond has a translation of (a,0,0) between the
% interacting spins. After assigning a matrix to this bond, the matrix of
% the second equivalent bond is generated using the (t,R) symmetry
% operator ( t - translation, R - rotation). The R rotation (90 degree
% roation around the *c*-axis) is applied to the J interactions matrix:

sq.addmatrix('label','J','value',[0 0 1;0 0 0;-1 0 0])
sq.addcoupling('J',1)
J = sq.couplingtable(1).matrix

%%
% In the above example a Dzyaloshinskii-Moriya interaction is defeind. The
% coupling matrix on the second equivalent bond J(:,:,2) is rotated by 90
% degree by the following operation: LATEXJ'=R\cdot J\cdot R^TPATEX. This can be also seen on the plot.

plot(sq,'aCoupling',false)







