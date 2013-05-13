function res = sw_test_sym()
% test symmetry on tetragonal structure
% 4-fold axis along z

tetra = sw;
tetra.genlattice('lat_const',[6 6 5],'sym','P 4');

tetra.addatom('r',[1/4 1/4 0]);

tetra.addmatrix('label',{'J1'});

tetra.gencoupling;
tetra.addcoupling('J1',1);

plot(tetra,'range',[-0.1 2.1; -0.1 2.1;-0.1 2.1]);

end