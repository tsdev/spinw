%%

model = sw_model('triAF',1);
plot(model)

% calculate LSWT spectrum
spec  = model.spinwave({[0 0 0] [1 1 0] 601});
EvBin = linspace(0,5,501);
EvCenter = (EvBin(2:end)+EvBin(1:(end-1)))/2;
spec  = sw_egrid(spec,'Evect',EvBin);

% FWHM of the Lorentzian component
lorFWHM = 0.2;
% should be larger than zero
gauFWHM = 1e-5;

spec = sw_instrument(spec,'dE',gauFWHM,'func',@(x,p)swfunc.voigtfwhm(x,[p lorFWHM 0]));

figure
subplot(2,1,1)
sw_plotspec(spec,'mode','color')
axis([0 1 0 5])
legend off

% plot cut at a selected Q value
Q = 1/2; % (1/2,1/2,0)

% find the index of the Q in the calculated Q points
idx = findmin(abs(spec.hkl(1,:)-Q));
% intensity
I = spec.swConv(:,idx);


subplot(2,1,2)
plot(EvCenter,I)
xlabel('Energy Transfer (meV)')
ylabel('Intensity (arb. units)')

title(sprintf('Q=(%4.2f,%4.2f,%4.2f)',[Q Q 0]))



