%% Ferromagnetic first neighbor Kagome lattice
% We define the kagome lattice using the "P -3" space group (Nr. 147),
% using the 3 fold symmetry, we need to define only one magnetic atom. In
% this space group all fisrt neighbor couplings will be equivalent (related
% by symmetry). Symmetry equivalent positions are automatically generated
% by the sw.atom() function. The magnetic atoms is Cu+ with S=1 spin.

FMkagome = spinw; 
FMkagome.genlattice('lat_const',[6 6 5],'angled',[90 90 120],'spgr','P -3')
FMkagome.addatom('r', [1/2 0 0], 'S', 1, 'label','MCu1','color','r')
display('Atomic positions as columns:') 
FMkagome.atom.r 
plot(FMkagome,'range',[2 2 1],'zoom',-0.5)

%% Create FM bonds
% The first neighbor bonds will be ferromagnetic, J = -1 meV. The
% sw.gencoupling() will use the space group operators to identify
% equivalent couplings, if two bonds have the same length but not symmetry
% related, then they will be identified as different bonds.

FMkagome.gencoupling('maxDistance',4) 
display('Rows: dlx, dly, dlz, at1, at2, idx, ma1, ma2, ma3')
FMkagome.couplingtable.table 
FMkagome.coupling 
display('Bond vectors (first three rows) and bond distances')
FMkagome.couplingtable.bondv 

FMkagome.addmatrix('label','J1','value',-1,'color','orange'); 
FMkagome.addcoupling('mat','J1','bond',1);
plot(FMkagome,'range',[2 2 1],'zoom',1.2)

%% FM magnetic structure
% All spins are paralle, pointing along the y-axis (perpendicular to ac
% plane). We use the "helical" mode of the sw.gencoupling() function, even
% though the structure is not helical. However in this mode all missing
% spins will be automatically created using the k-vector and normal vector
% and assuming helical magnetic structure. Thus if we give k = (0 0 0) and
% only the direction of the first spin in the unit cell, the code will
% create all other spin parallel to the first.

FMkagome.genmagstr('mode','helical','k',[0 0 0],'n',[0 1 0],'S',[0 1 0])
display('Magnetic structure with spins 1 2 ... as columns, xyz as rows:')
FMkagome.magstr
FMkagome.magstr.S
display('Ground state energy before optimization')
FMkagome.energy
plot(FMkagome,'range',[2 2 1])

%% Spin wave dispersion
% We calculate the spin wave dispersion. There are three modes, equal to
% the number of magnetic atom in the magnetic unit cell, denoted by
% different colors. At the zone center, the dispersion of the goldstone
% mode is parabolic, due to the FM interactions.

fmkSpec = FMkagome.spinwave({[-1/2 0 0] [0 0 0] [1/2 1/2 0] 100},'hermit',false);
fmkSpec = sw_neutron(fmkSpec);
fmkSpec = sw_egrid(fmkSpec, 'Evect',linspace(0,6.5,100),'component','Sperp'); 
sw_plotspec(fmkSpec,'mode',1,'colorbar',false,'axLim',[0 8])

%% Powder spectrum
% We plot the powder spectrum two different ways. First as calculated (to
% show the very strong Van Hoove singularity at the top of the dispersion),
% secondly convolute with a Gaussian along energy.

fmkPow = FMkagome.powspec(linspace(0,2.5,100),'Evect',linspace(0,7,250),'nRand',1000,'hermit',false);
figure;
sw_plotspec(fmkPow,'colorbar',true,'axLim',[0 0.05])
figure;
sw_plotspec(fmkPow,'colorbar',true,'axLim',[0 0.05],'dE',0.25,'norm',true)

%%
%  Written by
%  Bjorn Fak & Sandor Toth
%  06-June-2014