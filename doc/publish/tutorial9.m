%% k=0 Kagome antiferromagnet with DM interaction
% We create the lattice with 'P -3' space group and magnetic Cu+ with S=1
% spin.

DMkag = sw;
DMkag.fileid(0)
DMkag.genlattice('lat_const',[6 6 40],'angled',[90 90 120],'sym','P -3')
DMkag.addatom('r', [1/2 0 0],'S',1,'label', 'Cu1', 'color','r')
plot(DMkag,'range',[2 2 1],'zoom',-0.5)

%% Create bonds and Hamiltonian
% Generate the list of bonds and assign a Heisenberg exchange and weak
% Dzyaloshinskii-Moriya interaction to the fisrt neighbor bonds. The DM
% interaction vector can be easily created using the sw.setmatrix()
% function. The corresponding sw.getmatrix() function determines the
% symmetry allowed components of the DM vector. On the structure plot, the
% DM interaction vector is symbolized by a vector in the middle of the
% bond, pointing in the direction of the DM vector.

DMkag.gencoupling('maxDistance',7)
DMkag.addmatrix('label','DM1','value',1,'color','b')
DMkag.addmatrix('label','J1','value',1,'color','g')
DMkag.addcoupling('DM1',1)
DMkag.addcoupling('J1',1)
DMkag.setmatrix('label','DM1','pref',{[0 0 0.08]});
plot(DMkag,'range',[3 3 1],'zoom',-0.8)
display('J1='); DMkag.matrix.mat(:,:,2)
display('DM1='); DMkag.matrix.mat(:,:,1)

%% Generate magnetic structure
% We create a k = (0 0 0) magnetic structure, with the three spin directions
% in the unit cell (120 degree between neighbors). The spin vector
% components are given in the coordinate system of the lattice vectors
% (abc).

S0 = [1 -2 1; 2 -1 -1; 0 0 0];
DMkag.genmagstr('mode','direct','k',[0 0 0],'n',[0 0 1],'unitS','lu', 'S',S0); 
display('Ground state energy (meV/spin)')
DMkag.energy
plot(DMkag,'range',[3 3 1])

%% Calculate spin wave dispersion

Qlist = {[-1/2 0 0] [0 0 0] [1/2 1/2 0] 100};
dmkSpec = DMkag.spinwave(Qlist,'hermit',false);
sw_plotspec(dmkSpec,'mode',1,'axLim',[0 3],'colorbar',false,'colormap',[0 0 0],'dashed',true)


%% Powder spectrum
% The flat mode that is the zero energy mode lifted by the DM interaction
% is well visible on the powder spectrum.

dmkPow = DMkag.powspec(linspace(0,2.5,150),'Evect',linspace(0,3,250),'nRand',1000,'hermit',false);
figure;
sw_plotspec(dmkPow,'axLim',[0 0.5],'dE',0.02)

%%
%  Written by
%  Bjorn Fak & Sandor Toth
%  06-June-2014