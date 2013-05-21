function res = sw_test_sym()
% test symmetry on tetragonal structure
% 4-fold axis along z

%%
tetra = sw;
tetra.genlattice('lat_const',[6 6 5],'sym','P 4','angle',[90 90 90]*pi/180);

%tetra.addatom('r',[1/4 1/4 0]);
tetra.addatom('r',[0/2 0/2 0]);

tetra.addmatrix('label',{'J1'});
tetra.addmatrix('label',{'J2'},'color',[0; 255; 0]);
tetra.addmatrix('label',{'A'},'mat',diag([1 0 0]));

tetra.gencoupling;
tetra.addcoupling('J1',1);
tetra.addcoupling('J2',2);


tetra.addaniso('A',1);
%tetra.matrix.mat(:,:,3) = [0 1 0;1 0 0;0 0 0];
tetra.setmatrix('label','A','fid',1,'pref',{[1 1]})

%plot(tetra,'range',[2 2 1]);
end