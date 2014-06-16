%% Spin wave disperion of YVO3
% We compare our results with the model from: C. Ulrich, et al. PRL 91, 257202 (2003).
% see [[http://prl.aps.org/abstract/PRL/v91/i25/e257202]]
% We create crystal structure of YVO3 in the pseudocubic unit cell, doubled
% along c-axis. The magnetic atoms are V4+ with spin quantum number S=1/2.

a = 5.2821;
b = 5.6144;
c = 7.5283;

yvo3 = sw;
yvo3.fileid(0)
yvo3.genlattice('lat_const', [a/sqrt(2) b/sqrt(2) c])
yvo3.addatom('r',[0 0 0],'label','MV4','S',1/2,'color','gray')
yvo3.addatom('r',[0 0 1/2],'label','MV4','S',1/2,'color','gray')
yvo3.gencoupling
%yvo3.newcell({[1 1 0] [-1 1 0] [0 0 1]})
plot(yvo3)

%% Magnetic Hamiltonian
% The exchange constants are taken from the paper.

Jab   = 2.6;
Jc    = 3.1;
delta = 0.35;
K1    = 0.90;
K2    = 0.97;
d     = 1.15;

Jc1mat = -Jc*(1+delta)*eye(3) + diag([K2 0 0]) - [0 0 d;0 0 0;-d 0 0];
Jc2mat = -Jc*(1-delta)*eye(3) + diag([K2 0 0]) + [0 0 d;0 0 0;-d 0 0];

yvo3.addmatrix('label','Jab','value',Jab,'color','b')
yvo3.addmatrix('label','Jc1','value',Jc1mat,'color','YellowGreen')
yvo3.addmatrix('label','Jc2','value',Jc2mat,'color','purple')
yvo3.addmatrix('label','K1','value',-diag([K1 0 0]),'color','orange')
yvo3.addcoupling('Jab',[1 3])
yvo3.addcoupling('Jc1',2,2)
yvo3.addcoupling('Jc2',2,1)

yvo3.addaniso('K1')
plot(yvo3)

%theta = 1/2*atan(2*d/(2*Jc-K1-K2));

%% Magnetic structure
% We define G-type antiferromagnetic structure.

yvo3.genmagstr('mode','helical','S',[1 0 0],'nExt',[2 2 1],'k',[1/2 1/2 1],'n',[0 1 0])
plot(yvo3,'ellAniso',1.5)

%%
Optimising magnetic structure, assuming it is planar, the meaning of the x parameters are the following: (phi1, phi2, ... phiN, kx, ky kz, nTheta, nPhi):
par_opt.xmin = [zeros(1,8)    , 0 0 0, 0  0 ];
par_opt.xmax = [ones(1,8)*2*pi, 0 0 0, pi pi];
par_opt.func = @gm_planar;
par_opt.nRun = 10;
optRes = yvo3.optmagstr(par_opt);
Check canting angles (compare phi to theta):
M   = yvo3.mag_str.S;
phi = atan2(M(1,:),M(3,:))*180/pi;
Calculating spin wave dispersion:
specYVO3 = yvo3.spinwave({[3/4 3/4 0] [1/2 1/2 0] [1/2 1/2 1] });
specYVO3 = sw_neutron(specYVO3);
specYVO3 = sw_conv(specYVO3,'Evect',linspace(0,25,nE));
Plotting of spin wave spectra with imaginary components:
figure
subplot(3,1,1);
sw_plotspec(specYVO3,'mode',1,'aHandle',gca,'imag',true);
subplot(3,1,2);
sw_plotspec(specYVO3,'mode',2,'aHandle',gca,'imag',true);
subplot(3,1,3);
sw_plotspec(specYVO3,'mode',3,'aHandle',gca);

%%
%  Written by
%  Sandor Toth
%  16-June-2014