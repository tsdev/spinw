%% Optimize chain in magnetic field

chain = spinw;
chain.genlattice('lat_const',[3 4 4])
chain.addatom('r',[0 0 0],'S',1)
chain.addmatrix('label','J1','value',1)
chain.addmatrix('label','D','value',[0 0 0.1])
chain.gencoupling
chain.addcoupling('mat','J1','bond',1)
chain.addcoupling('mat','D','bond',1)
plot(chain)


%% Optimize chain in magnetic field

chain = spinw;
chain.genlattice('lat_const',[6 4 4])
chain.addatom('r',[0 0 0],'S',1)
chain.addatom('r',[1/2 0 0],'S',1)
chain.addmatrix('label','J1','value',1)
chain.addmatrix('label','D','value',[0 0 0.1])
chain.gencoupling
chain.addcoupling('mat','J1','bond',1)
chain.addcoupling('mat','D','bond',1)
plot(chain)

%% k-vector

res = chain.optmagk;
chain.energy
