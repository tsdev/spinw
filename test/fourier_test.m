%% kagome lattice

kag = spinw;
kag.genlattice('lat_const',[6 6 4],'angled',[90 90 120])
kag.addatom('r',[1/2 0   0],'S',1)
kag.addatom('r',[1/2 1/2 0],'S',1)
kag.addatom('r',[0   1/2 0],'S',1)

%kag.newcell('bvect',{[2 0 0] [0 1 0] [0 0 1]})
kag.quickham(1)
%plot(kag)

%% calculate fourier transformation

for ii = 1:10:100
    Q = reshape(sw_qgrid('bin',{[0 0.01 2] [0 0.01 2] [1 ii]}),3,[]);
    tic;
    F2 = kag.fourier(Q,'isomode','auto');
    t1(ii) = toc;
end

