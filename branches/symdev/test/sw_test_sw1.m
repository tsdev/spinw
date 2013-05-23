%% define crystal structure of Na2IrO3

% data from: Choi, S. K., Coldea, et al. (2012). PRL, 108(12), 127204
% space group C2/m, crystallographic parameters at 300 K

nairo = sw;
%nairo.genlattice('lat_const',[5.427 9.395 5.614],'angle',[90 109.037 90]*pi/180,'sym',12);
nairo.genlattice('lat_const',[5.427 9.395 5.614],'angle',[90 109.037 90]*pi/180,'sym','C 2/m');

% add magnetic Ir
nairo.addatom('r',[1/2; 0.167; 0],'S',1/2,'label',{'Ir'},'color',[140; 28; 22]);
nairo.addatom('r',[0 1/2 1/2;0 0 0.340; 0 1/2 1/2],'S',[0 0 0],'label',{'Na1' 'Na2' 'Na3'},'color',ones(3)*233);
nairo.addatom('r',[0.748 0.711; 0.178 0; 0.789 0.204],'S',[0 0],'label',{'O1', 'O2'},'color',[45 45;118 118;125 125]);

plot(nairo)

%% define prototype magnetic Hamiltonian

% regenerate crystal with P1 symmetry
nairo.newcell({[1 0 0] [0 1 0] [0 0 1]})

% Kitaev terms
JKxx=zeros(3); JKxx(1,1)=1;
JKyy=zeros(3); JKyy(2,2)=1;
JKzz=zeros(3); JKzz(3,3)=1;

% rotation needs to be added, see sw_rot
nairo.addmatrix('mat',JKxx,'color',[255; 0; 0],'label','JKxx');
nairo.addmatrix('mat',JKyy,'color',[0; 255; 0],'label','JKyy');
nairo.addmatrix('mat',JKzz,'color',[0; 0; 255],'label','JKzz');

% Heisenberg terms
nairo.addmatrix('mat',eye(3),'color',[128; 128; 128],'label','J1');
nairo.addmatrix('mat',eye(3),'color',[200; 128; 128],'label','J2');
nairo.addmatrix('mat',eye(3),'color',[128; 200; 128],'label','J3');

% generate couplings up to 8 Angstrom
nairo.gencoupling('maxdistance',8);

% add J1, J2 and J3 and JK couplings
nairo.addcoupling('J1',[1 2]);
nairo.addcoupling('J2',[3 4]);
nairo.addcoupling('J3',[7 8]);
% anisotropic Kitaev couplings
nairo.addcoupling('JKxx',1,[2 4]);
nairo.addcoupling('JKyy',1,[1 3]);
nairo.addcoupling('JKzz',2);

plot(nairo,'range',[-0.1 2.1;-0.1 2.1;-0.1 0.1])

%% define Q scans

nQ = 200;
nE = 400;
dQ = 0.764e-6;
Qp{1} = [ -1;   0; 0]+dQ;
Qp{2} = [  0;   0; 0]+dQ;
Qp{3} = [  0;   1; 0]+dQ;
Qp{4} = [  1;   1; 0]+dQ;
Qp{5} = [1/2; 1/2; 0]+dQ;
Qp{6} = [  0;   0; 0]+dQ;

%% case d stripy order
% energy per spin is -0.2912

J1 =  1;
J2 =  0;
J3 =  0;
JK =  1.33;
Na2IrO3fun(nairo,[J1 J2 J3 JK]);

opt_par.xmin = [0 0 0 0,        0 0 0, 0 0];
opt_par.xmax = [[1 1 1 1]*2*pi, 0 0 0, pi 0];
opt_par.func = @gm_planar;
opt_par.nRun = 10;

optRes = nairo.optmagstr(opt_par);

% calculate spin waves with zig-zag order
E  = linspace(0,4,nE);

specD = nairo.spinwave([Qp {nQ}]);

specD = sw_neutron(specD,'pol',false);
specD = sw_conv(specD,'convmode','Sxx','evect',E);

figure
sw_plotspec(specD,'mode',3,'ahandle',gca,'imag',false,'convE',0.05,'axLim',0.5);
sw_plotspec(specD,'mode',1,'ahandle',gca,'imag',false,'dashed',true,'colorbar',false);


%% case e-f stripy order
% energy per spin is -0.375

J1 =  1;
J2 =  0;
J3 =  0;
JK =  2;
Na2IrO3fun(nairo,[J1 J2 J3 JK]);

opt_par.xmin = [0 0 0 0,        0 0 0, 0 0];
opt_par.xmax = [[1 1 1 1]*2*pi, 0 0 0, pi 0];
opt_par.func = @gm_planar;

nTry = 10;
Emin = 0;

for ii = 1:nTry
    optRes = nairo.optmagstr(opt_par);
    if optRes.e < Emin
        Emin = optRes.e;
        xmin = optRes.x;
    end
end

nairo.genmagstr('mode','func','func',@gm_planar,'x0',xmin);
plot(nairo,'range',[-0.1 1.1;-0.1 1.1;-0.1 0.1],'scaleS',2)

% calculate spin waves with zig-zag order
E  = linspace(0,4,nE);

specEF = nairo.spinwave([Qp {nQ}],'evect',E,'polfact',false);

figure
sw_plotspec(specEF,'mode',3,'ahandle',gca,'imag',false,'convE',0.05);
sw_plotspec(specEF,'mode',1,'ahandle',gca,'imag',false,'nocolmap',true,'dashed',true);


%% case g stripy order
% optimise magnetic structure assuming planar
% --> finds tripy order

J1 =  1;
J2 =  0.26;
J3 = -0.2;
JK = 0;
Na2IrO3fun(nairo,[J1 J2 J3 JK]);

opt_par.xmin = [0 0 0 0,        0 0 0, 0 0];
opt_par.xmax = [[0 1 1 1]*2*pi, 0 0 0, 0 0];
opt_par.func = @gm_planar;

nairo.optmagstr(opt_par);

% calculate spin waves with stripy order
E  = linspace(0,2,nE);

specG = nairo.spinwave([Qp {nQ}],'evect',E);

figure
sw_plotspec(specG,'mode',3,'ahandle',gca,'imag',false,'convE',0.05);
sw_plotspec(specG,'mode',1,'ahandle',gca,'imag',false,'nocolmap',true,'dashed',true);

%% case h zig-zag order
% optimise magnetic structure assuming planar
% --> finds tripy order

J1 =  1;
J2 =  0.78;
J3 =  0.9;
JK =  0;
Na2IrO3fun(nairo,[J1 J2 J3 JK]);

opt_par.xmin = [0 0 0 0,        0 0 0, 0 0];
opt_par.xmax = [[0 1 1 1]*2*pi, 0 0 0, 0 0];
opt_par.func = @gm_planar;

nairo.optmagstr(opt_par);

% calculate spin waves with zig-zag order
E  = linspace(0,4,nE);

specH = nairo.spinwave([Qp {nQ}],'evect',E);

figure
sw_plotspec(specH,'mode',3,'ahandle',gca,'imag',false,'convE',0.05);
sw_plotspec(specH,'mode',1,'ahandle',gca,'imag',false,'nocolmap',true,'dashed',true);

%% powder spectra
J1 =  4.17;
J2 =  0.78*J1;
J3 =  0.9*J1;
JK =  0;
Na2IrO3fun(nairo,[J1 J2 J3 JK]);

hklA = linspace(0.2,1.5,30);

specHpow = nairo.powspec(hklA,'evect',linspace(0,8,100),'nRand',1e4);

sw_plotspec(specHpow,'convE',0.5)

%% case ij zig-zag order
% optimise magnetic structure assuming planar

% DISPERSION DEPENDS ON THE SPIN DIRECTION
% THE M POINT IS NOT UNIQUE

J1 =  1;
J2 =  0.23;
J3 =  0.51;
JK =  1.33;
Na2IrO3fun(nairo,[J1 J2 J3 JK]);

opt_par.xmin = [0 0 0 0,        0 0 0, 0 0];
opt_par.xmax = [[0 1 1 1]*2*pi, 0 0 0, pi 0];
opt_par.func = @gm_planar;

nairo.optmagstr(opt_par);

% calculate spin waves with zig-zag order
E  = linspace(0,3,nE);

specIJ = nairo.spinwave([Qp {nQ}]);
specIJ = sw_neutron(specIJ,'pol',false);

convmode = {'Sxx' 'Szz'};
figure
for ii = 1:2
    subplot(2,1,ii);
    specIJ = sw_conv(specIJ,'convmode',convmode{ii},'evect',E);
    
    sw_plotspec(specIJ,'mode',3,'ahandle',gca,'imag',false,'convE',0.05);
    sw_plotspec(specIJ,'mode',1,'ahandle',gca,'imag',false,'dashed',true,'colorbar',false,'title',false);
    caxis([0 0.5])
end










