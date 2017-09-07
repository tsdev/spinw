%% alpha-CuV2O7

cuvo = spinw('alpha-Cu2V2O7.cif');
cuvo.addmatrix('label','J1','value',2.67);
cuvo.addmatrix('label','J2','value',2.99);
cuvo.addmatrix('label','J3','value',5.42);
cuvo.addmatrix('label','G','value',diag([1 -1 -1])*0.282);
cuvo.addmatrix('label','D','value',[2.79 0 0]);

cuvo.gencoupling
cuvo.addcoupling('bond',1,'mat','D')
cuvo.addcoupling('bond',1,'mat','G')
cuvo.addcoupling('bond',1,'mat','J1')
cuvo.addcoupling('bond',2,'mat','J2')
cuvo.addcoupling('bond',3,'mat','J3')
plot(cuvo,'atommode','mag','range',[1 2 2])

% k=0 magnetic structure
cuvo.optmagsteep('nrun',1e4)
cuvo.energy

%% spinw wave

spec = cuvo.spinwave({[0 1 0] [0 3 0] 501});
spec = sw_egrid(spec);
figure
sw_plotspec(spec)
legend off

%% low field phase

figure

cuvo.field([6 0 0])
cuvo.optmagsteep('nRun',1e4)

spec = cuvo.spinwave({[0 1 0] [0 3 0] 501});
spec = sw_egrid(spec,'Evect',linspace(0,4.5,501));
subplot(2,2,1)
sw_plotspec(spec,'dE',0.5)
axis([0 2 0 4.5])
legend off
title('B = +6 T')
cuvo.field([-6 0 0])
cuvo.optmagsteep('nRun',1e4)

spec = cuvo.spinwave({[0 1 0] [0 3 0] 501});
spec = sw_egrid(spec,'Evect',linspace(0,4.5,501));
subplot(2,2,2)
sw_plotspec(spec,'dE',0.5)
axis([0 2 0 4.5])
legend off
title('B = -6 T')

%% high field phase
cuvo.field([0 0 0])
cuvo.addmatrix('label','G','value',0);

% get incommensurate propagation vector
cuvo.optmagk
% generate random magnetic structure
cuvo.genmagstr('mode','helical','k',[0 cuvo.magstr.k(2) 0],'S',[zeros(1,16);rand(2,16)],'n',[1 0 0])

res = cuvo.optmagsteep('nRun',1e4)
cuvo.addmatrix('label','G','value',diag([1 -1 -1])*0.282);
cuvo.energy

spec = cuvo.spinwave({[0 1 0] [0 3 0] 501});
spec = sw_egrid(spec,'Evect',linspace(0,4.5,501));
subplot(2,2,3)
sw_plotspec(spec,'dE',0.5)
axis([0 2 0 4.5])
legend off
title('B = +10 T')
cuvo.field([-6 0 0])
cuvo.optmagsteep('nRun',1e4)










