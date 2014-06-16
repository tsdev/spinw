% Spin wave simulation for LuMnO3
% parameters from PRL 111, 257202 (2013)
% "Magnon Breakdown in a Two Dimensional Triangular Lattice Heisenberg Antiferromagnet of Multiferroic LuMnO3"
%
%
%
%% crystal structure

importciffile = true;

if importciffile
    % import crystal structure from .cif file, .cif file can be found in
    % ICSD Karlsruhe database
    lumno = sw('LuMnO3.cif');
else
    % imput crystal structure directly
    lumno = sw;
    lumno.genlattice('lat_const',[6.038 6.038 11.361],'angled',[90 90 120],'sym','P 63 c m')
    lumno.addatom('r',[-0.3355 -0.3355 -0.00077],'label','Mn3+','S',2)
end

%% define couplings

seti = 2;

switch seti
    case 1
        % set 1 of couplings
        J1 = 9; J2 = 1.4; J3 = 0; J1c = 0.018; J2c = 0; D1 = 0.28; D2 = 0.006;
    case 2
        % set 2 of couplings
        J1 = 6.4; J2 = 1.3; J3 = -0.15; J1c = -0.009; J2c = 0.009; D1 = 0.5; D2 = -0.009;
    otherwise
        error('Wrong set!');
end

% generate list of bonds based on symmetry
lumno.gencoupling

% define coupling and anisotropy matrices
lumno.addmatrix('label','J1' ,'color',[255 0   0],  'value',J1)
lumno.addmatrix('label','J2' ,'color',[0   255 0],  'value',J2)
lumno.addmatrix('label','J3' ,'color',[0   0   255],'value',J3)
lumno.addmatrix('label','J1c','color',[255 180 0],  'value',J1c)
lumno.addmatrix('label','J2c','color',[0   180 255],'value',J2c)
lumno.addmatrix('label','D1','color',[255 0   180],'value',diag([0 0 D1]))

% assign matrices to bonds
lumno.addcoupling('J1' ,2)
lumno.addcoupling('J2' ,1)
lumno.addcoupling('J3' ,6)
lumno.addcoupling('J1c',4)
lumno.addcoupling('J2c',3)

% assign D1 anisotropy to all magnetic atoms
lumno.addaniso('D1')

% D2 is not included, since the direction of this additional anisotropy is
% undefined in the paper

% plot only magnetic atoms
plot(lumno,'pNonMagAtom',false,'range',[1 1 1])

%% magnetic structure

%optimise the magnetic structure based on classical energy
% assuming planar structure in the ab plane
optPar.func = @gm_planar;
% limits of the optimised parameters, see also optRes.xname list
%              Phi1 Phi2 Phi3 Phi4 Phi5 Phi6 kx ky kz nTheta nPhi
optPar.xmin = [0    0    0    0    0    0    0  0  0  0      0];
optPar.xmax = [0    2*pi 2*pi 2*pi 2*pi 2*pi 0  0  0  0      0];
% number of runs with random initial parameters
optPar.nRun = 10;

% optimisation
optRes = lumno.optmagstr(optPar);

% check the Phi angles and give the exact values below
phiOpt = optRes.x(1:6)*180/pi;

% exact values are simply rounded from the appriximate ones to degree
phi = round(phiOpt)*pi/180;

% save magnetic structure into lumno object with the exact values
lumno.genmagstr('mode','func','func',@gm_planar,'x0',[phi 0 0 0 0 0])


%% calculate spin waves

lSpec = lumno.spinwave({[-3/4 -1/4 0] [-1/2 -1/2 0] [0 0 0] [1 0 0] [1/2 1/2 0] [1/2 0 0] [0 1 0] [0 1 1]},'hermit',false);

% calculate neutron scattering cross section and save it inot lSpec.Sperp
lSpec = sw_neutron(lSpec);
% create energy grid from the neutron scattering intensities for plotting
lSpec = sw_egrid(lSpec,'component','Sperp');

%% plot spectrum

% plot dispersion only, 12 modes is 2x6, where the 6 invisible ones are in
% negative energy
figure
subplot(2,1,1)
sw_plotspec(lSpec,'mode',1,'linestyle','o','dashed',true,'colormap',@jet,'dE',0.8,'axlim',[0 30],'aHandle',gca,'imag',true)

% plot color spectrum
subplot(2,1,2)
sw_plotspec(lSpec,'mode',3,'linestyle','-','dashed',true,'colormap',@fireprint,'dE',0.8,'axlim',[0 30],'aHandle',gca)








