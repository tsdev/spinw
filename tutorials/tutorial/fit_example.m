%% Generate artificial data
% We use a single gaussian function with linear background and simulating
% Poisson statistics.

% parameters
Amp  = 70;
Cen  = 50;
Wid  = 20;
Bkg  = 20;

% ideal function values
xDat = 0:100;
pDat = [Amp Cen Wid 0 Bkg];
yDat = fitfun.ngaussbkg(xDat,pDat);

% store in sig
sig = struct;
sig.x = xDat;
% generate Poisson variables
for ii = 1:numel(xDat)
    sig.y(ii) = poisson(yDat(ii));
    sig.e(ii) = sqrt(sig.y(ii));
end

% fit data with slightly bad starting values and calculate the error of the
% simulated function at 95% confidence level.
[pFit,yCalc,fitRes3] = ndbase.lm(sig,@fitfun.ngaussbkg,[300 40 10 0 15],'confLev',0.95,'extraStat',true);

figure
subplot(2,2,1)
% Plot the gridded parameter search as a function of the first 2 parameter
% values.

Ngrid   = 51;
[p1,p2] = ndgrid(linspace(0,150,Ngrid),linspace(25,75,Ngrid));
Chi2    = zeros(Ngrid);

for ii = 1:Ngrid
    for jj = 1:Ngrid
        yTemp = fitfun.ngaussbkg(sig.x,[p1(ii,jj) p2(ii,jj) Wid 0 Bkg]);
        dy = sig.y - yTemp;
        Chi2(ii,jj) = dy*dy'/(numel(sig.x)-numel(pFit));
    end
end

% plot mesh
mesh(p1,p2,log10(Chi2))
xlabel('p_1')
ylabel('p_2')
zlabel('log_{10}(\chi^2_\nu)','rotation',90)

subplot(4,2,2)
% Plot the parameter values during the fitting process.
maxIt = fitRes3.cvgHst(end,1);
% number of parameters
Np = numel(pDat);
% plot the parameter values
plot(fitRes3.cvgHst(:,1),fitRes3.cvgHst(:,(2:Np)+1), '-o','linewidth',4);
for ii=1:Np
    text(1.02*fitRes3.cvgHst(end,1),fitRes3.cvgHst(end,1+ii), sprintf('%d',ii) );
end
xlim([1 maxIt])
ylabel('p values')

subplot(4,2,4)
% plot the chi2 and lambda values
semilogy(fitRes3.cvgHst(:,1) , [fitRes3.cvgHst(:,Np+2) fitRes3.cvgHst(:,Np+3)], '-o','linewidth',4)
text(0.8*fitRes3.cvgHst(end,1),0.2*fitRes3.cvgHst(end,Np+2), '\chi^2_\nu','FontSize',16 );%,'color','blue'
text(0.8*fitRes3.cvgHst(end,1),2.5*fitRes3.cvgHst(end,Np+3), '\lambda','FontSize',16 );%,'color',dkgrn
ylabel('\chi^2_\nu, \lambda')
xlabel('function calls')
ylim([1e-5 1e3])
xlim([1 maxIt])

subplot(2,2,3)
% Plot te data and the fit function and the 95% confidence levels.
errorbar(sig.x,sig.y,sig.e)
hold on
plot(sig.x,fitfun.ngaussbkg(sig.x,fitRes3.p),'-r','linewidth',2)
hold on
plot(sig.x,fitfun.ngaussbkg(sig.x,fitRes3.p)+fitRes3.sigY,'.k','markersize',8)
plot(sig.x,fitfun.ngaussbkg(sig.x,fitRes3.p)-fitRes3.sigY,'.k','markersize',8)

subplot(2,2,4)
% Plot histogram of yDat-yCalc, it should be gaussian
histogram(sig.y-yCalc,10)
xlabel('y-f(x,p)')
ylabel('#')