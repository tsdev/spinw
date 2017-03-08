%% antiferromagnetic chain

afc = spinw;
afc.genlattice('lat_const',[3 4 4])
afc.addatom('r',[  0 0 0],'S',1)
afc.addmatrix('label','A','value',diag([0 0 0.1]))
afc.addmatrix('label','J1','value',1)
afc.addmatrix('label','J2-','value',1/3,'color','orange')
afc.gencoupling
afc.addcoupling('mat','J1', 'bond',1)
afc.addcoupling('mat','J2-','bond',5)
afc.addaniso('A')
%afc.optmagsteep;
plot(afc,'range',[2 1 1])

%% opt magnetic structure

%afc.field([0 0 7])
afc.optmagstr('func',@gm_spherical3d,'xmin',[0 0, 0 0 0,0 0],'xmax',[pi/2 0,1/2 0 0,0 0])
E0 = afc.energy;
%% magnetic field

optRes = afc.optmagsteep('random',false,'nRun',400);
%plot(afc)
figure;
plot(optRes.e,'o-')

%% spin wave
spec = afc.spinwave({[0 0 0] [2 0 0] 500},'hermit',true);
spec = sw_egrid(spec,'Evect',linspace(0,3,300));
figure
sw_plotspec(spec);%,'mode','disp','imag',true);
axis([0 2 0 3])


