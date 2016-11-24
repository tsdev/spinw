%% antiferromagnetic chain

afc = sw;
afc.genlattice('lat_const',[4 4 6])
afc.addatom('r',[  0 0 0],'S',1)
afc.addatom('r',[1/2 0 0],'S',1)
afc.addmatrix('label','A','value',diag([0 0 -0.1]))
afc.addmatrix('label','J','value',1)
afc.gencoupling
afc.addcoupling('J',1)
afc.addaniso('A')
afc.optmagsteep;
plot(afc)

%% magnetic field

afc.field([0 0 7])
afc.optmagsteep('random',true,'nRun',400);
plot(afc)

%% spin wave
spec = afc.spinwave({[0 0 0] [2 0 0] 500},'hermit',true);
spec = sw_egrid(spec,'Evect',linspace(0,3,300));
figure
sw_plotspec(spec);%,'mode','disp','imag',true);
axis([0 2 0 3])


