%% generate structure with 1st neghbor bonds

ca2ruo4 = spinw('/Users/sandortoth/Documents/structures/Ca2RuO4/ca2ruo4.cif');
ca2ruo4.unit_cell.S(2) = 1;
plot(ca2ruo4)

ca2ruo4.gencoupling()
ca2ruo4.addmatrix('label','J','value',1)
ca2ruo4.addcoupling('bond',1,'mat','J')

%% symmetry allowed elements of the exchange interaction

ca2ruo4.getmatrix('mat','J')

%% assign values to the exchange matrix

J = 5.2;
a = 0.10*J;
e = 1;
E = 21.5;
A = 1.0;

Jmat = [J A 0;A J 0;0 0 J-a];
Amat = [e 0 0; 0 0 0; 0 0 E];
ca2ruo4.addmatrix('label','J','value',Jmat)
ca2ruo4.addmatrix('label','A','value',Amat)
ca2ruo4.addaniso('A')
% show bond directional interactions, generated using the symmetry
% operators of the space group
%plot(ca2ruo4,'range',[2 2 1])

% optimize magnetic structure

ca2ruo4.genmagstr('mode','random','S',[1 0 0]','k',[0 0 0])
magRes = ca2ruo4.optmagsteep('nRun',1e4);
ca2ruo4.energy
% calculate spin wave dispersion

A = [1/2 0 0];
M = [1/2 1/2 0];
B = [1 0 0];
G = [0 0 0];
spec = ca2ruo4.spinwave({A M B G M});
%spec = ca2ruo4.spinwave({[1/2 1/2 0] [1 0 0] [1 1 0] [0 0 0] [1 0 0]});
spec = sw_egrid(spec,'component','Sxx+Syy+Szz','Evect',linspace(0,65,501));
spec = sw_instrument(spec,'dE',4.2);
figure
sw_plotspec(spec,'qlabel',{'A' 'M' 'B' '\Gamma' 'M'},'axLim',[0 65],'mode','disp','colormap',[0 0 0],'dashed',false,'linestyle','--')
hold on
sw_plotspec(spec,'qlabel',{'(1/2,0)' '(1/2,1/2)' '(1,0)' '(0,0)' '(1/2,1/2)'},'axLim',[0 0.05],'mode','color')
legend('off')
colorbar
title('')


