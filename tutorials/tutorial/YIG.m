%% Model YIG spin wave spectrum
% We will reproduce the calculation of the YIG spin wave spectrum to
% compare to this paper:
% PRL 117, 217201 (2016)    http://journals.aps.org/prl/abstract/10.1103/PhysRevLett.117.217201
% arXiv:1607.03263          https://arxiv.org/abs/1607.03263

% get the .cif file from Google Drive, the file is also available on ICSD
yig = spinw('https://goo.gl/2ez9Ys');

% color differently the two Fe sublattice
yig.unit_cell.color(:,3) = swplot.color('r');
yig.unit_cell.color(:,2) = swplot.color('b');

% spin quantum number of Fe3+ ions, determined automatically by SpinW
S0 = max(yig.unit_cell.S);

% normalize spins to S=1 as it is in the paper
yig.unit_cell.S = yig.unit_cell.S/S0;

hFig = plot(yig,'atomMode','mag');
swplot.plotchem('atom1','Fe1','atom2','O','limit',6)
swplot.plotchem('atom1','Fe2','atom2','O','limit',4,'replace',false)
 

%% show primitive cell

plot(yig,'atomMode','mag')

% new basis vectors in rows
pBV = [1/2 1/2 -1/2;-1/2 1/2 1/2;1/2 -1/2 1/2];
% lattice constant of YIG
lat = yig.abc(1);

swplot.plot('type','arrow','position',cat(3,zeros(3),pBV'))

%% generate the bonds
% An interesting symmetry property of YIG in the "I a -3 d" space group is
% that the bonds type 3 and type 4 have the exact same length, however they
% are not related by symmetry. This can be easily seen by checking the
% center psition of the bonds:
yig.gencoupling('maxDistance',6);
yig.getmatrix('bond',3);
yig.getmatrix('bond',4);

% The (7/8,1/8,1/8) position belongs to the 48g Wyckoff position, while the
% (7/8,7/8,7/8) position is 16b. Thus the exchange interactions on the two
% bonds can be different, even though previous models of YIG assumed they
% are equal.


%% create spin Hamiltonian
% change from BCC to primitive cubic cell
T = yig.newcell('bvect',pBV);

% exchange values from the paper
Jad = sw_converter(9.60e-21,'J','THz','photon');
Jdd = sw_converter(3.24e-21,'J','THz','photon');
Jaa = sw_converter(0.92e-21,'J','THz','photon');

% scale the interactions from classical moment size to quantum model
Scl = sqrt(S0*(S0+1));
yig.quickham([Jad Jdd Jaa]/Scl)

% add external field and convert from the standard SpinW unit (meV) to THz
yig.field([0 0 0.01]*sw_converter(1,'meV','THz','photon'))
yig.optmagsteep
yig.genmagstr('mode','rotate','n',[0 0 1])
yig.optmagsteep

if sum(yig.mag_str.F(3,:),2)<0
    yig.mag_str.F = -yig.mag_str.F;
end

plot(yig,'atomMode','mag','atomLegend',false)

%% Spin wave dispersion to compare with the paper

Q0  = [1 2 3]';
Q_N = T*([ 1/2  1/2    0]'+Q0);
Q_G = T*([   0    0    0]'+Q0);
Q_H = T*([   0    0    1]'+Q0);

spec = yig.spinwave({Q_N Q_G Q_H 501});

%% plot the spin wave dispersion

spec = sw_egrid(spec,'component','Sxy-Syx','Evect',linspace(0,28,501));
spec = sw_instrument(spec,'dE',0.75);

figure
sw_plotspec(spec,'mode','disp','colormap',[0 0 0])
hold on
sw_plotspec(spec,'mode','color','imag',true)
colormap(sw_cbrewer('RdBu'))
title('YIG low temperature spin wave spectrum')
ylabel('Energy (THz)')

legend off
colorbar off
set(gca,'XTickLabel',{'N' '\Gamma' 'H'})
caxis([-0.05 0.05])

%% powder spactrum
pref = swpref;
Q = linspace(0,4,101);
E = linspace(0,30,501);
nRand = 1e2;
pref.tid = 2;

spec = yig.powspec(Q,'Evect',E,'nRand',nRand);

%%
tri  = sw_model('triAF',1);
Q = linspace(0,4,101);
E = linspace(0,30,501);
nRand = 1e3;
pref.tid = 0;
profile on
spec = tri.powspec(Q,'Evect',E,'nRand',nRand);
profile off

