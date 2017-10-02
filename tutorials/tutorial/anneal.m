% Test to reproduce results of this paper:
% P. Lacorre, et al., J. of Magn. and Magn. Mat., 71, 63(1987).

%% Defines lattice for antiferromagnetic chains

% AFM chain along a-axis
hChain = sw_model('chain',1);

% Adds D easy-axis single-ion anisotropy
hChain.addmatrix('value',diag([0 0 -0.2]),'label','A','color','red')

% Adds the defined D anisotropy to all magnetic atoms
hChain.addaniso('A')

% Removes units
hChain.unit.kB = 1;
hChain.unit.muB = 1;

%% Plots crystal structure

plot(hChain,'range',[1 1/2 1/2]);

%% Random initial spin configuration with 150x1x1 unit cell

hChain.genmagstr('mode','random','nExt',[150 1 1])

%% Plots initial magnetic structure

plot(hChain,'range',[10 1/2 1/2]);

%% Perform simulated annealing with periodic boundary condition

anneal_param.boundary  = {'per' 'free' 'free'};
anneal_param.verbosity = 2;
anneal_param.cool      = @(T)(0.8*T);
anneal_param.initT     = 40;
anneal_param.endT      = 1e-3;
anneal_param.nMC       = 1000;
aRes = hChain.anneal(anneal_param);

%% Plots final magnetic structure

plot(hChain,'range',[4 1/2 1/2]);

%% Simulated annealing with periodic boundary conditions in magnetic field

hChain.field([0 0 0]);

param = struct('verbosity',1,'cool',@(T)(0.8*T),'initT',40,'endT',1e-2,'nMC',1000,'nStat',1,'random',true);
hChain.anneal(param);

fieldSweep = [linspace(0,5,10) linspace(5,0,10)];
loopp = struct('x',fieldSweep,'func',@(obj,x)obj.field([0 0 x]),'nMC',2e3,'nStat',1e3,'verbosity',1);

profile on
pStat = hChain.annealloop(loop_param);
profile off

M = squeeze(sum(pStat.avgM,2))/hChain.nmagext;

%% Plot M versus external magnetic field

figure;
plot(fieldSweep,M(3,:),'r.-')

title(sprintf('S=1 Heisenberg antiferromagnetic chain, 150x1x1 unit cell, T=%g K',param.endT));
xlabel('Magnetic Field (T)')
ylabel('Magnetic moment')
legend('M_z')
axis([0 5 0 1.2])