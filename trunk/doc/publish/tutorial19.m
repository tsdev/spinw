%% Cu and Fe chain

chain = sw;
chain.genlattice('lat_const',[3 8 4],'sym','P 1')
chain.addatom('label','MCu2','r',[0 0 0])
chain.addatom('label','MFe2','r',[0 1/2 0])

chain.gencoupling
chain.addmatrix('label','JCu','value',1,'color','r')
chain.addmatrix('label','JFe','value',1,'color','b')
chain.addmatrix('label','JCuFe','value',-0.1,'color','orange')

chain.addcoupling('JCu',1)
chain.addcoupling('JFe',2)
chain.addcoupling('JCuFe',[5 6])
plot(chain)

%% magnetic structure

chain.genmagstr('mode','direct','S',[0 0;1 1;0 0],'k',[1/2 0 0])
chain.genmagstr('mode','helical','nExt',[2 1 1])
plot(chain)

%% correlation function with form factor

spec = chain.spinwave({[0 0 0] [1 0 0] 500},'formfact',true);
spec = sw_egrid(spec,'component','Sxx+Syy+Szz');
figure
sw_plotspec(spec,'mode','color','dE',0.2)
figure
sw_plotspec(spec,'mode','disp','imag',true)

%% only copper correlations

spec = chain.spinwave({[0 0 0] [1 0 0] 500},'formfact',{1 0});
spec = sw_egrid(spec,'component','Sxx+Syy+Szz');
figure
sw_plotspec(spec,'mode','color','dE',0.2)

%% Fe intensity

FMchain = sw_model('chain',-1);
FMchain.genmagstr('mode','helical','nExt',[1 1 1])
%FMchain.addtwin('axis',[0 0 1],'phid',90)

spec = FMchain.spinwave({[0 0 0] [20 0 0] 101},'formfact',{'MIr4'});
spec = sw_egrid(spec,'component','Sxx+Syy+Szz','sumtwin',false);
figure
F = sw_mff('MIr4',spec.hklA);

%hold all
%plot(linspace(0,12,300),F.^2)
plot(sqrt(sum(spec.hklA.^2,1)),sum(spec.swInt(1:end/2,:),1)./F.^2)
%plot(sqrt(sum(spec.hklA.^2,1)),sum(spec.swInt(1:end/2,:),1))

