%% Spin-spin correlation function of two different coupled magnetic ion
% In this tutorial, we demonstrate how the spin wave spectrum of two
% different coupled ion can be calculated. Using the magnetic form factor,
% the spin-spin correlation function on the two ions can be separated. This
% method is also usefull for simulating RIXS (Resonant Inelastic X-ray
% Scattering) on spin waves, that typically measures only the spin waves on
% the resonant ion.

%% Cu and Fe chain
% We define a unit cell with an Fe2+ and Cu2+ ions, the size of the spin is
% determined by the code using the information stored in the  ion.dat file.
% We couple antiferromagnetically the Fe and Cu atoms in separate chains
% running along the a-axis. The chains are weakly coupled along the b-axis.

FeCuChain = sw;
FeCuChain.fileid(0)
FeCuChain.genlattice('lat_const',[3 8 4],'sym','P 1')
FeCuChain.addatom('label','MCu2','r',[0 0 0])
FeCuChain.addatom('label','MFe2','r',[0 1/2 0])

FeCuChain.gencoupling
FeCuChain.addmatrix('label','J_{Cu-Cu}','value',1,'color','r')
FeCuChain.addmatrix('label','J_{Fe-Fe}','value',1,'color','b')
FeCuChain.addmatrix('label','J_{Cu-Fe}','value',-0.1,'color','orange')

FeCuChain.addcoupling('J_{Cu-Cu}',1)
FeCuChain.addcoupling('J_{Fe-Fe}',2)
FeCuChain.addcoupling('J_{Cu-Fe}',[5 6])
plot(FeCuChain,'range',[2 0.7 1],'zoom',1)

%% Define magnetic structure
% The classical ground state of the above define Hamiltonian is
% antiferromagnetic chains ordered ferromagnetically along the b-direction.
% Here we define a magnetic unit cell that is 2x1x1 of the crystal unit
% cell to incorporate the antiferromagnetic structure.

FeCuChain.genmagstr('mode','direct','S',[0 0;1 1;0 0],'k',[1/2 0 0])
FeCuChain.genmagstr('mode','helical','nExt',[2 1 1])
plot(FeCuChain,'range',[2 0.7 1],'zoom',1,'rSpin',0.08,'lHeadSpin',0.75,'sSpin',1.5,'zoom',2)

%% Spin-spin correlation function including the form factor
% We calculate the spin wave intensity (withouth the polarisation factor)
% and include the magnetic form factor of both the Cu2+ and Fe2+. The spin
% wave dispersion plot includes the imaginary part of the spin wave
% energies as well (plotted as red dots) to check that the magnetic ground
% state is right.

spec = FeCuChain.spinwave({[0 0 0] [1 0 0] 501},'formfact',true,'hermit',false);
spec = sw_egrid(spec,'component','Sxx+Syy+Szz');
figure
sw_plotspec(spec,'mode','disp','imag',true,'axLim',[0 4.5])
figure
sw_plotspec(spec,'mode','color','dE',0.2)

%% Copper only correlations
% To see only the spin-spin correlations on the Cu2+ ions, the
% sw.spinwave() function has an option 'formfact', this option allows to
% change separately the form factor function of each magnetic ion in the
% unit cell separately. Since we have two magnetic ions, {1 0} means that
% the Cu2+ (defined first) will have a form factor of identity, while the
% Fe2+ ion will have form factor of zero (zero intensity). As the resulting
% plot shows, the upper branch comes from mostly the Fe2+ chain. This is
% resonable, since the spin wave energy is ~J*S and the Fe2+ ions have
% larger S. There is also a small mixing around the zone center and
% boundaries due to the weak J_{Cu-Fe} coupling.

spec = FeCuChain.spinwave({[0 0 0] [1 0 0] 503},'formfact',{1 0},'hermit',false);
spec = sw_egrid(spec,'component','Sxx+Syy+Szz');
figure
sw_plotspec(spec,'mode','color','dE',0.2)

%% Iron only correlations
% By changing the form factor option value to {0 1} the spin-spin
% correlations on the Fe ions are plotted. It is clearly visible, that
% the goldstone mode at (1/2,0,0) has a strong Fe2+ component as well.

spec = FeCuChain.spinwave({[0 0 0] [1 0 0] 500},'formfact',{0 1});
spec = sw_egrid(spec,'component','Sxx+Syy+Szz');
figure
sw_plotspec(spec,'mode','color','dE',0.2)

%% Reduced moments
% The size of the reduced moment due to the zero point energy of the spin
% waves can be calculated as well. The calculated momen reduction is
% finite, that shows that within the spin wave approach, the 2D coupled
% plane of Heisenberg spins at zero temperature is ordered. The moment
% reduction is much larger on the spin-2 Fe ion (dS = 0.577) than on the
% spin-1/2 (dS = 0.137) Cu ions. This is due to the fact, that for the Fe
% ions the interchain coupling is much weaker than the chain coupling, this
% leads to large moment reduction, while the Cu ions the interchain
% coupling is much larger relative to the intrachain exchange.

FeCuChain.moment('nRand',1e4)

%%
%  Written by
%  Sandor Toth
%  2-September-2014

