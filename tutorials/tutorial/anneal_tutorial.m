% Test to reproduce results of this paper:
% P. Lacorre, et al., J. of Magn. and Magn. Mat., 71, 63(1987).

%% Defines lattice for antiferromagnetic chains

% AFM chain along a-axis
hChain = sw_model('chain',2);

% Adds D easy-axis single-ion anisotropy
hChain.addmatrix('value',diag([0 0 -0.2]),'label','A','color','red')

% Adds the defined D anisotropy to all magnetic atoms
hChain.addaniso('A')

% Removes units
hChain.unit.kB = 1;
hChain.unit.muB = 1;

%% Plots crystal structure

plot(hChain,'range',[1 1/2 1/2]);


%% Plots initial magnetic structure

plot(hChain,'range',[10 1/2 1/2]);

%% Perform simulated annealing with periodic boundary condition in magnetic field
% 150 unit cell along a-axis

T = 1e-2;
hChain.field([0 0 0]);
param = struct('verbosity',2,'cool',@(T)(0.8*T),'initT',40,'endT',T,'nMC',1e3,'nStat',0,'random',true,'nExt',[150 1 1]);
aRes = hChain.anneal(param);

fieldSweep = [linspace(0,5,40) linspace(5,0,40)];
loop_param = struct('x',fieldSweep,'func',@(obj,x)obj.field([0 0 x]),'nMC',2e4,'nStat',1e4,'tid',2,'fineT',T);

pStat = hChain.annealloop(loop_param);


M = squeeze(sum(pStat.avgM,2))/hChain.nmagext;

%% Plot M versus external magnetic field

hold on
plot(pStat.x,M(3,:),'g.-')

title(sprintf('S=1 Heisenberg antiferromagnetic chain, 150x1x1 unit cell, T=%g K',param.endT));
xlabel('Magnetic Field (T)')
ylabel('Magnetic moment')
legend('M_z')
axis([0 5 0 1.2])