%% Spin wave spectrum of La2CuO4
% Crystal structure of La2CuO4 contains Cu2+ atoms with S = 1/2 spin.
% [1] R. Coldea, Phys. Rev. Lett. 86, 5377 (2001).

J   = 138.3;
Jp  = 2;
Jpp = 2;
Jc  = 38;

lacuo = sw_model('squareAF',[J-Jc/2 Jp-Jc/4 Jpp]);
lacuo.unit_cell.S = 1/2;

%%

Qlist = {[3/4 1/4 0] [1/2 1/2 0] [1/2 0 0] [3/4 1/4 0] [1 0 0] [1/2 0 0] 100};
lacuoSpec = lacuo.spinwave(Qlist,'hermit',false);
lacuoSpec = sw_neutron(lacuoSpec);
lacuoSpec = sw_egrid(lacuoSpec);
figure
sw_plotspec(lacuoSpec,'mode',1,'axLim',[0 600],'dashed',true)
figure
lacuoSpec = sw_omegasum(lacuoSpec,'zeroint',1e-5,'tol',1e-3);
sw_plotspec(lacuoSpec,'mode',2,'axLim',[0 20],'dashed',true,'colormap',[0 0 0])


%% Spin wave spectrum of La2CuO4
% Crystal structure of La2CuO4 contains Cu2+ atoms with S = 1/2 spin.

lacuo = sw;
sw_addsym('x+1/2,y,z+1/2;-x+1/2,y+1/2,-z;x+1/2,-y+1/2,z;-x,-y,-z','B m a b');
lacuo.genlattice('lat_const',[5.3561 5.3888 13.173],'sym','B m a b')
lacuo.addatom('label','MCu2','r',[0 0 0])
lacuo.addatom('label','O','r',[1/4 1/4 -0.0055])
lacuo.addatom('label','O','r',[0 0.0267 0.1832])
lacuo.addatom('label','La','r',[0 -0.0055 0.3609])
plot(lacuo)

%% Hamiltonian

J   = 138.3;
Jp  = 2;
Jpp = 2;
Jc  = 38;

lacuo.addmatrix('label','J',  'value',J-Jc/2, 'color','r')
lacuo.addmatrix('label','Jp', 'value',Jp-Jc/4,'color','g')
lacuo.addmatrix('label','Jpp','value',Jpp,    'color','b')

lacuo.gencoupling
lacuo.addcoupling('J',  1)
lacuo.addcoupling('Jp', 2)
lacuo.addcoupling('Jpp',3)

%% Magnetic structure

x1 = [0 0 0 0, 0 0 0, 0 0];
x2 = [0 [1 1 1]*2*pi, 1 1 0, 0 0];
optRes = lacuo.optmagstr('func',@gm_planar,'xmin',x1,'xmax',x2,'nRun',10);

Qlist = {[3/4 1/4 0] [1/2 1/2 0] [1/2 0 0] [3/4 1/4 0] [1 0 0] [1/2 0 0] 100};
lacuoSpec = lacuo.spinwave(Qlist,'hermit',false);
sw_plotspec(lacuoSpec,'mode',1)











