%%

tri1 = sw_model('triAF',1);
tri2 = tri1.copy;



tic
spec1 = tri1.spinwave({[0 0 0] [1 1 0] 5001});
spec1 = sw_egrid(spec1);
toc

tic
spec2 = tri2.spinwavefast({[0 0 0] [1 1 0] 5001});
spec2 = sw_egrid(spec2);
toc

figure
subplot(2,1,1)
sw_plotspec(spec1,'mode','color','dE',0.2)
subplot(2,1,2)
sw_plotspec(spec2,'mode','color','dE',0.2)