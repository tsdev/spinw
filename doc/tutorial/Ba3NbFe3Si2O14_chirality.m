% Ba3NbFe3Si2O14 chiral compound
% crystal structure: PRL 101, 247201 (2008)
% spin wave:         PRL 106, 207201 (2011)

%% crystal structure

banb = sw;
banb.genlattice('lat_const',[8.539 8.539 5.2414],'angle',[90 90 120]*pi/180,'sym','P 3 2 1');

banb.addatom('r',[0.24964;0;1/2],'S',5/2,'color',[128;128;128],'label','MFe3');

banb.gencoupling;

J1 = 0.85;
J2 = 0.24;
J3 = 0.053;
J4 = 0.017;
J5 = 0.24;

banb.addmatrix('mat',eye(3)*J1,'label','J1','color',[150;150;255])
banb.addmatrix('mat',eye(3)*J2,'label','J2','color',[0;0;0])
banb.addmatrix('mat',eye(3)*J3,'label','J3','color',[0;200;0])
banb.addmatrix('mat',eye(3)*J4,'label','J4','color',[0;0;255])
banb.addmatrix('mat',eye(3)*J5,'label','J5','color',[125;0;125])

banb.addcoupling('J1',1)
banb.addcoupling('J4',2)
banb.addcoupling('J2',3)
banb.addcoupling('J5',4)
banb.addcoupling('J3',5)

plot(banb,'range',[2 2 2])

%% determine magnetic structure
% minimum energy -5.4522 meV/spin

banb.genmagstr('mode','direct','S',[1 -1/2 -1/2;0 -sqrt(3)/2 sqrt(3)/2;0 0 0],'next',[1 1 1])
xmin = [0 0 0        0 0 0 0 0];
xmax = [0 2*pi*[1 1] 0 0 1 0 0];

optRes = banb.optmagstr('func',@gm_planar,'xmin',xmin,'xmax',xmax,'nRun',10);

Emin = banb.energy;

%% fix magnstic structure 
% using the chirality values for triangle and helix

eD = -1;
eH = +1;
banb.genmagstr('mode','helical','S',[1 -1/2 -1/2;0 eD*sqrt(3)/2 -eD*sqrt(3)/2;0 0 0],'k',[0 0 eH*1/7],'n',[0 0 1],'next',[1 1 7])


%% plot magnetic structure

% eD
plot(banb,'range',[-0.5 0.5;-0.5 0.5;-0.1 1.1])

% eH
plot(banb,'range',[0.1 0.9; 0.1 0.9; 0 7],'coplanar',1)

%% spin waves -- epsilon_T = -1 crystal structure

eH = +1;
eD = -1;
banb.genmagstr('mode','helical','S',[1 -1/2 -1/2;0 eD*sqrt(3)/2 -eD*sqrt(3)/2;0 0 0],'k',[0 0 eH*1/7],'n',[0 0 1],'next',[1 1 7])

banbSpec = banb.spinwave({[0 -1 1] [0 -1 -2] 250});
banbSpec.hkl = -banbSpec.hkl;
banbSpec.hklA = -banbSpec.hklA;
banbSpec = sw_neutron(banbSpec,'pol',true,'n',[1 0 0]);
banbSpec = sw_egrid(banbSpec,'component','Sperp','Evect',linspace(0,6,500));

figure
sw_plotspec(banbSpec,'mode',3,'dE',0.25,'ahandle',gca);
caxis([0 1])
title('\epsilon_T = -1','fontsize',16)
colormap(jet)

banbSpec = banb.spinwave({[0 -1 1] [0 -1 -2] 300});
banbSpec.hkl  = -banbSpec.hkl;
banbSpec.hklA = -banbSpec.hklA;
banbSpec = sw_neutron(banbSpec,'pol',true,'n',[1 0 0]);
banbSpec = sw_egrid(banbSpec,'component','Myz-Mzy','Evect',linspace(0,6,500));

figure
subplot(2,1,1)
sw_plotspec(banbSpec,'mode',3,'dE',0.25,'ahandle',gca,'imag',true);
caxis([-1 1])
colormap(makecolormap([0 1 0],[1 1 1],[0 0 1],1000));
title('\epsilon_T = -1 (\epsilon_H = +1, \epsilon_\Delta = -1)','fontsize',16)

eH = -1;
eD = +1;
banb.genmagstr('mode','helical','S',[1 -1/2 -1/2;0 eD*sqrt(3)/2 -eD*sqrt(3)/2;0 0 0],'k',[0 0 eH*1/7],'n',[0 0 1],'next',[1 1 7])

banbSpec = banb.spinwave({[0 -1 1] [0 -1 -2] 300});
banbSpec.hkl  = -banbSpec.hkl;
banbSpec.hklA = -banbSpec.hklA;
banbSpec = sw_neutron(banbSpec,'pol',true,'n',[1 0 0]);
banbSpec = sw_egrid(banbSpec,'component','Myz-Mzy','Evect',linspace(0,6,500));

subplot(2,1,2)
sw_plotspec(banbSpec,'mode',3,'dE',0.25,'ahandle',gca,'imag',true);
caxis([-1 1])
colormap(makecolormap([0 1 0],[1 1 1],[0 0 1],1000));
title('\epsilon_T = -1 (\epsilon_H = -1, \epsilon_\Delta = +1)','fontsize',16)


%% change crystal chirality to epsilon_T = +1

% change the chirality of the crystal
banb.setmatrix('label','J3','pref',{J5})
banb.setmatrix('label','J5','pref',{J3})

%% spin waves -- epsilon_T = -1 crystal structure

eH = +1;
eD = +1;
banb.genmagstr('mode','helical','S',[1 -1/2 -1/2;0 eD*sqrt(3)/2 -eD*sqrt(3)/2;0 0 0],'k',[0 0 eH*1/7],'n',[0 0 1],'next',[1 1 7])

banbSpec = banb.spinwave({[0 -1 1] [0 -1 -2] 250});
banbSpec.hkl = -banbSpec.hkl;
banbSpec.hklA = -banbSpec.hklA;
banbSpec = sw_neutron(banbSpec,'pol',true,'n',[1 0 0]);
banbSpec = sw_egrid(banbSpec,'component','Sperp','Evect',linspace(0,6,500));

figure
sw_plotspec(banbSpec,'mode',3,'dE',0.25,'ahandle',gca);
caxis([0 1])
title('\epsilon_T = -1','fontsize',16)
colormap(jet)

banbSpec = banb.spinwave({[0 -1 1] [0 -1 -2] 300});
banbSpec.hkl  = -banbSpec.hkl;
banbSpec.hklA = -banbSpec.hklA;
banbSpec = sw_neutron(banbSpec,'pol',true,'n',[1 0 0]);
banbSpec = sw_egrid(banbSpec,'component','Myz-Mzy','Evect',linspace(0,6,500));

figure
subplot(2,1,1)
sw_plotspec(banbSpec,'mode',3,'dE',0.25,'ahandle',gca,'imag',true);
caxis([-1 1])
colormap(makecolormap([0 1 0],[1 1 1],[0 0 1],1000));
title('\epsilon_T = -1 (\epsilon_H = +1, \epsilon_\Delta = -1)','fontsize',16)

eH = -1;
eD = -1;
banb.genmagstr('mode','helical','S',[1 -1/2 -1/2;0 eD*sqrt(3)/2 -eD*sqrt(3)/2;0 0 0],'k',[0 0 eH*1/7],'n',[0 0 1],'next',[1 1 7])

banbSpec = banb.spinwave({[0 -1 1] [0 -1 -2] 300});
banbSpec.hkl  = -banbSpec.hkl;
banbSpec.hklA = -banbSpec.hklA;
banbSpec = sw_neutron(banbSpec,'pol',true,'n',[1 0 0]);
banbSpec = sw_egrid(banbSpec,'component','Myz-Mzy','Evect',linspace(0,6,500));

subplot(2,1,2)
sw_plotspec(banbSpec,'mode',3,'dE',0.25,'ahandle',gca,'imag',true);
caxis([-1 1])
colormap(makecolormap([0 1 0],[1 1 1],[0 0 1],1000));
title('\epsilon_T = -1 (\epsilon_H = -1, \epsilon_\Delta = +1)','fontsize',16)


%% constant energy surface using Horace

horace_on;

%profile -memory on
d3dobj = d3d(banb.abc,[0 1 0 0],[-3.5,0.01,0.5],[0 0 1 0],[-2.5,0.01,0.5],[0 0 0 1],[0.5,0.025,6]);
d3dobj = disp2sqw_eval(d3dobj,@banb.horace,[],0.1);
%profile off
plot(d3dobj);




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% incommensurate
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% change crystal chirality to epsilon_T = -1

% change the chirality of the crystal
banb.setmatrix('label','J3','pref',{J3})
banb.setmatrix('label','J5','pref',{J5})








%% spin waves -- epsilon_T = -1 crystal structure

eH = +1;
eD = -1;
banb.genmagstr('mode','helical','S',[1 -1/2 -1/2;0 eD*sqrt(3)/2 -eD*sqrt(3)/2;0 0 0],'k',[0 0 eH*1/7],'n',[0 0 1],'next',[1 1 1])

banbSpec = banb.spinwave({[0 -1 1] [0 -1 -2] 250});
banbSpec.hkl = -banbSpec.hkl;
banbSpec.hklA = -banbSpec.hklA;
banbSpec = sw_neutron(banbSpec,'pol',true,'n',[1 0 0]);
banbSpec = sw_egrid(banbSpec,'component','Sperp','Evect',linspace(0,6,500));

figure
sw_plotspec(banbSpec,'mode',3,'dE',0.25,'ahandle',gca);
caxis([0 1])
title('\epsilon_T = -1','fontsize',16)
colormap(jet)

banbSpec = banb.spinwave({[0 -1 1] [0 -1 -2] 300});
banbSpec.hkl  = -banbSpec.hkl;
banbSpec.hklA = -banbSpec.hklA;
banbSpec = sw_neutron(banbSpec,'pol',true,'n',[1 0 0]);
banbSpec = sw_egrid(banbSpec,'component','Myz-Mzy','Evect',linspace(0,6,500));

figure
subplot(2,1,1)
sw_plotspec(banbSpec,'mode',3,'dE',0.25,'ahandle',gca,'imag',true);
caxis([-1 1])
colormap(makecolormap([0 1 0],[1 1 1],[0 0 1],1000));
title('\epsilon_T = -1 (\epsilon_H = +1, \epsilon_\Delta = -1)','fontsize',16)

eH = -1;
eD = +1;
banb.genmagstr('mode','helical','S',[1 -1/2 -1/2;0 eD*sqrt(3)/2 -eD*sqrt(3)/2;0 0 0],'k',[0 0 eH*1/7],'n',[0 0 1],'next',[1 1 1])

banbSpec = banb.spinwave({[0 -1 1] [0 -1 -2] 300});
banbSpec.hkl  = -banbSpec.hkl;
banbSpec.hklA = -banbSpec.hklA;
banbSpec = sw_neutron(banbSpec,'pol',true,'n',[1 0 0]);
banbSpec = sw_egrid(banbSpec,'component','Myz-Mzy','Evect',linspace(0,6,500));

subplot(2,1,2)
sw_plotspec(banbSpec,'mode',3,'dE',0.25,'ahandle',gca,'imag',true);
caxis([-1 1])
colormap(makecolormap([0 1 0],[1 1 1],[0 0 1],1000));
title('\epsilon_T = -1 (\epsilon_H = -1, \epsilon_\Delta = +1)','fontsize',16)


%% change crystal chirality to epsilon_T = +1

% change the chirality of the crystal
banb.setmatrix('label','J3','pref',{J5})
banb.setmatrix('label','J5','pref',{J3})

%% spin waves -- epsilon_T = -1 crystal structure

eH = +1;
eD = +1;
banb.genmagstr('mode','helical','S',[1 -1/2 -1/2;0 eD*sqrt(3)/2 -eD*sqrt(3)/2;0 0 0],'k',[0 0 eH*1/7],'n',[0 0 1],'next',[1 1 1])

banbSpec = banb.spinwave({[0 -1 1] [0 -1 -2] 250});
banbSpec.hkl = -banbSpec.hkl;
banbSpec.hklA = -banbSpec.hklA;
banbSpec = sw_neutron(banbSpec,'pol',true,'n',[1 0 0]);
banbSpec = sw_egrid(banbSpec,'component','Sperp','Evect',linspace(0,6,500));

figure
sw_plotspec(banbSpec,'mode',3,'dE',0.25,'ahandle',gca);
caxis([0 1])
title('\epsilon_T = -1','fontsize',16)
colormap(jet)

banbSpec = banb.spinwave({[0 -1 1] [0 -1 -2] 300});
banbSpec.hkl  = -banbSpec.hkl;
banbSpec.hklA = -banbSpec.hklA;
banbSpec = sw_neutron(banbSpec,'pol',true,'n',[1 0 0]);
banbSpec = sw_egrid(banbSpec,'component','Myz-Mzy','Evect',linspace(0,6,500));

figure
subplot(2,1,1)
sw_plotspec(banbSpec,'mode',3,'dE',0.25,'ahandle',gca,'imag',true);
caxis([-1 1])
colormap(makecolormap([0 1 0],[1 1 1],[0 0 1],1000));
title('\epsilon_T = -1 (\epsilon_H = +1, \epsilon_\Delta = -1)','fontsize',16)

eH = -1;
eD = -1;
banb.genmagstr('mode','helical','S',[1 -1/2 -1/2;0 eD*sqrt(3)/2 -eD*sqrt(3)/2;0 0 0],'k',[0 0 eH*1/7],'n',[0 0 1],'next',[1 1 1])

banbSpec = banb.spinwave({[0 -1 1] [0 -1 -2] 300});
banbSpec.hkl  = -banbSpec.hkl;
banbSpec.hklA = -banbSpec.hklA;
banbSpec = sw_neutron(banbSpec,'pol',true,'n',[1 0 0]);
banbSpec = sw_egrid(banbSpec,'component','Myz-Mzy','Evect',linspace(0,6,500));

subplot(2,1,2)
sw_plotspec(banbSpec,'mode',3,'dE',0.25,'ahandle',gca,'imag',true);
caxis([-1 1])
colormap(makecolormap([0 1 0],[1 1 1],[0 0 1],1000));
title('\epsilon_T = -1 (\epsilon_H = -1, \epsilon_\Delta = +1)','fontsize',16)



