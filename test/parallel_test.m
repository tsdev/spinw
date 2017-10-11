%%

tri1 = sw_model('triAF',1);
tri2 = tri1.copy;
tri2.addtwin('axis',[0 0 1],'phid',[60 120])
tri2.genmagstr('mode','helical','k',[1/3 1/3 0],'n',[0 0 1],'S',[1 0 0]','nExt',[3 3 1]);

%Q = {[0 0 0] [1 1 0] 501};
Q = [0.5 0.5 0]';

spec1 = tri1.spinwave(Q);
spec1 = sw_egrid(spec1,'Evect',linspace(0,4,501));

spec2 = spinwavefast(tri2,Q);
spec2 = sw_egrid(spec2,'Evect',linspace(0,4,501));

%%
figure
subplot(2,1,1)
sw_plotspec(spec1,'mode','color','dE',0.2)
subplot(2,1,2)
sw_plotspec(spec2,'mode','color','dE',0.2)