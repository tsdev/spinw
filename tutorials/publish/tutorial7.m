%% k=0 Kagome antiferromagnet
% We create a lattice with space group "P -3" where all first neighbor
% bonds are symmetry equivalent and add a magnetic Cr+ with S=1 spin.

AFkagome = spinw;
AFkagome.fileid(0)
AFkagome.genlattice('lat_const',[6 6 10],'angled',[90 90 120],'spgr','P -3')
AFkagome.addatom('r',[1/2 0 0],'S', 1,'label','MCu1','color','r')
plot(AFkagome,'range',[2 2 1],'zoom',-0.5)

%% Create bonds
% Generate the list of bonds and lists them.

AFkagome.gencoupling('maxDistance',7) 
display('Rows: dlx, dly, dlz, at1, at2, idx, ma1, ma2, ma3')
nBond = size(AFkagome.couplingtable([1 2]).table,2);
fprintf(['' repmat('%4d ',1,nBond) '\n'],AFkagome.couplingtable([1 2]).table');
display('Bond vectors (first three rows) and bond distances')
AFkagome.couplingtable([1 2]).bondv 

%% Hamiltonian
% We create AFM first neighbor interaction and weak 2nd neighbor AFM
% exchange.

AFkagome.addmatrix('label','J1','value',1.00,'color','r')
AFkagome.addmatrix('label','J2','value',0.11,'color','g')
AFkagome.addcoupling('mat','J1','bond',1)
AFkagome.addcoupling('mat','J2','bond',2)
plot(AFkagome,'range',[3 3 1],'zoom',-0.5)

%% Generate magnetic structure
% We create a k = (0 0 0) magnetic structure, with the three spin directions
% in the unit cell (120 degree between neighbors). The spin vector
% components are given in the coordinate system of the lattice vectors
% (abc).

S0 = [1 -2 1; 2 -1 -1; 0 0 0];
AFkagome.genmagstr('mode','direct','k',[0 0 0],'n',[0 0 1],'unitS','lu','S',S0); 
display('Magnetic structure with spins 1 2 ... as columns, xyz as rows:')
AFkagome.magstr
AFkagome.magstr.S
display('Magnetic atoms as columns:')
AFkagome.magtable.R
display('Magnetic spins:')
AFkagome.magtable.M
display('Ground state energy (meV/spin)')
AFkagome.energy
plot(AFkagome,'range',[3 3 1])

%% Calculate spin wave dispersion

afkSpec = AFkagome.spinwave({[-1/2 0 0] [0 0 0] [1/2 1/2 0] 100},'hermit',false);
sw_plotspec(afkSpec,'mode',1,'axLim',[0 3],'colorbar',false,'colormap',[0 0 0],'dashed',true)

%% Powder spectrum

afkPow = AFkagome.powspec(linspace(0,2.5,150),'Evect',linspace(0,3,250),'nRand',1000,'hermit',false);
figure;
sw_plotspec(afkPow,'axLim',[0 0.2])

%%
%  Written by
%  Bjorn Fak & Sandor Toth
%  06-June-2014