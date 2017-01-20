%% Ferromagnetic kagome lattice 
% We create the kagome lattice with up to 4th neighbor interactions. The
% symmetry related atoms are denoted by MCu1(i)_j, where i is the index of
% independent atomic positions, j is the index of the generated atomic
% positions of the i-th independent position.

kagome4 = spinw; 
kagome4.fileid(0)
kagome4.genlattice('lat_const',[6 6 40],'angled',[90 90 120],'spgr','P -3');
kagome4.addatom('r', [1/2 0 0],'S', 1,'label','MCu1','color','r');
display('Atomic positions as columns:') 
kagome4.atom.r
plot(kagome4)

%% Define Hamiltonian
% We add couplings up to 4th neighbor interactions. If the generation of
% the bond tables would depend on distance, J3 and Jd would be equivalent.
% However using the 'P -3' space group the two type of bonds are
% inequivalent, as physically expected in real systems (J3 goes through an
% intermediate magnetic atom, while Jd is not).

kagome4.gencoupling('maxDistance',7); 
display('Rows: dlx, dly, dlz, at1, at2, idx, ma1, ma2, ma3')
kagome4.couplingtable.table 
kagome4.coupling 
display('Bond vectors (first three rows) and bond distances')
kagome4.couplingtable.bondv 

kagome4.addmatrix('label','J1','value',-1.00,'color','g')
kagome4.addmatrix('label','J2','value', 0.10,'color','r')
kagome4.addmatrix('label','J3-','value', 0.00,'color','orange')
kagome4.addmatrix('label','Jd','value', 0.17,'color','b')

kagome4.addcoupling('mat','J1','bond',1);
kagome4.addcoupling('mat','J2','bond',2); 
kagome4.addcoupling('mat','J3-','bond',3); 
kagome4.addcoupling('mat','Jd','bond',4);
plot(kagome4,'range',[2 2 1],'zoom',-0.8)

%% Magnetic structure
% For strong FM 1str neighbour and weak further neighbor interaction the
% ground state is ferromagnetic.

kagome4.genmagstr('mode','helical','k',[0 0 0],'n',[0 1 0],'S',[0 1 0]); 
display('Magnetic structure with spins 1 2 ... as columns, xyz as rows:')
kagome4.mag_str
kagome4.mag_str.S
display('Ground state energy (meV/spin)')
kagome4.energy
plot(kagome4)

%% Spin wave dispersion

kag4Spec = kagome4.spinwave({[-1/2 0 0] [0 0 0] [1/2 1/2 0] 200},'hermit',false);
kag4Spec = sw_neutron(kag4Spec);
kag4Spec = sw_egrid(kag4Spec,'Evect',linspace(0,6.5,100)) ;
sw_plotspec(kag4Spec,'mode',1,'axLim',[0 8],'colorbar',false,'colormap',[0 0 0])

%% Powder averaged spectrum

kag4Pow = kagome4.powspec(linspace(0,2.5,100),'Evect',linspace(0,7,250),'nRand',1000,'hermit',false);
figure;
sw_plotspec(kag4Pow,'axLim',[0 0.2],'colorbar',true)

%%
%  Written by
%  Bjorn Fak & Sandor Toth
%  06-June-2014