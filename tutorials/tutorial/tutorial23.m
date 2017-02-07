%% Define crystal structure of Sr3Fe2O7

SrF = spinw;
SrF.genlattice('lat_const', [3.8 3.8 20.4],'sym','I 4/m m m')
SrF.addatom('r',[0 0 0.0972],'S',2,'label','MFe4','color','blue')
SrF.gencoupling

%% define magnetic Hamiltonian
% Value of the exchange constants

J1  = -7.20;
Jc1 = -5.10;
J2  =  1.05;
J3  =  2.10;
Jc3 =  0.01;
D   = diag([0 0 -0.06]);

SrF.addmatrix('label','J1', 'value',J1, 'color','red')
SrF.addmatrix('label','Jc1','value',Jc1,'color','orange')
SrF.addmatrix('label','J2', 'value',J2, 'color','white')
SrF.addmatrix('label','J3', 'value',J3, 'color','gray')
SrF.addmatrix('label','Jc3','value',Jc3,'color','yellow')
SrF.addmatrix('label','D',  'value',D,  'color','purple')

% to add the two ion coupling
SrF.addcoupling('mat','J1','bond',1)
SrF.addcoupling('mat','Jc1','bond',2)
SrF.addcoupling('mat','J2','bond',3)
SrF.addcoupling('mat','J3','bond',7)
SrF.addcoupling('mat','Jc3','bond',6)
SrF.addaniso('D')

% Define magnetic structure
S0 = [0 -1.162 0 -1.162; 0  1.162 0  1.162;2 -1.140 2 -1.140];
SrF.genmagstr('mode','helical','S',S0,'k',[1/7 1/7 1],'n',[1 1 0])
%SrF.optmagsteep
SrF.plot

%% Plot figure

figure
nHkl = 501;
Evect = 0:0.5:35;

% Fig. 3(f)
subplot(1,5,1)
specFerr = SrF.spinwave({[0.6 0.14 0] [1.4 0.14 0] nHkl},'hermit',false);
specFerr = sw_magdomain(specFerr,'axis',[0 0 1],'angled',[0 90 180 270]);
specFerr = sw_neutron(specFerr);
specFerr = sw_egrid(specFerr,'Evect',Evect);
specFerr = sw_instrument(specFerr,'Ei', 40,'dE',1.72,'dQ',0.05,'norm',false);
sw_plotspec(specFerr,'mode','color','legend',false,'axLim',[0 0.6]);
colormap(jet)

% Fig. 3(g)
subplot(1,5,2)
specFerr = SrF.spinwave({[0.6 -0.4 0] [1.4 0.4 0] nHkl},'hermit',false);
specFerr = sw_magdomain(specFerr,'axis',[0 0 1],'angled',[0 90 180 270]);
specFerr = sw_neutron(specFerr);
specFerr = sw_egrid(specFerr,'Evect',Evect);
specFerr = sw_instrument(specFerr,'Ei', 40,'dE',1.72,'dQ',0.05,'norm',false);
sw_plotspec(specFerr,'mode','color','axLim',[0 0.2],'legend',false)
colormap(jet)

% Fig. 3 (h)
subplot(1,5,3)
specFerr = SrF.spinwave({[-0.4 -0.4 5] [0.4 0.4 5] nHkl},'hermit',false);
specFerr = sw_magdomain(specFerr,'axis',[0 0 1],'angled',[0 90 180 270]);
specFerr = sw_neutron(specFerr);
specFerr = sw_egrid(specFerr,'Evect',Evect);
specFerr = sw_instrument(specFerr,'Ei',20,'dE',0.012,'dQ',0.0167,'norm',false);
sw_plotspec(specFerr,'mode',3,'axLim',[0 0.5],'legend',false)
colormap(jet)

% Fig. 3 (i)
subplot(1,5,4)
specFerr = SrF.spinwave({[-0.4 -0.4 5] [0.4 0.4 5] nHkl},'hermit',false);
specFerr = sw_magdomain(specFerr,'axis',[0 0 1],'angled',[0 90 180 270]);
specFerr = sw_neutron(specFerr);
specFerr = sw_egrid(specFerr,'Evect',Evect);
resMat = [0 10 20 30 40 50 60;6 4.5 3.25 2.12 1.15 0.4 0]';
specFerr = sw_instrument(specFerr,'Ei', 60,'dE',resMat,'dQ',0.05, 'norm',false);
sw_plotspec(specFerr,'mode',3,'axLim',[0 0.2],'legend',false)
colormap(jet)
axis([0.0 0.8 7 28]);

% Fig. 3 (j)
subplot(1,5,5)
specFerr = SrF.spinwave({[-0.4 -0.4 7] [0.4 0.4 7] nHkl},'hermit',false);
specFerr = sw_magdomain(specFerr,'axis',[0 0 1],'angled',[0 90 180 270]);
specFerr = sw_neutron(specFerr);
specFerr = sw_egrid(specFerr,'Evect',Evect);
specFerr = sw_instrument(specFerr,'Ei', 60,'dE',resMat,'dQ',0.05, 'norm',false);
sw_plotspec(specFerr,'mode','color','axLim',[0 0.2],'legend',false)
colormap(jet)
axis([0.0 0.8 8 35]);