%% simulate non-magnetic inpurities on the square lattice

% probability that the site is magnetic
p = 0.95;

S = 1/2;

N1 = 10;
N2 = 10;
nMagSite = rand(N1,N2)<p;

sq = sw;
sq.genlattice('lat_const',[N1*3 N2*3 5])

for ii = 1:N1
    for jj = 1:N2
        sq.addatom('r',[(ii-1)/N1 (jj-1)/N2 0],'S',nMagSite(ii,jj)*S)
    end
end

sq.addmatrix('label','J1','value',1)
sq.addmatrix('label','J2','value',0.1,'color','green')


sq.gencoupling
sq.addcoupling('J1',1)
sq.addcoupling('J2',2)
%sq.genmagstr('mode','random')
sq.plot('labelatom',false,'sspin',4)

%% monte carlo
% --> the ground state is trivial, since there is no frustration
% E1 = -0.3169

pMC.initT = 10;
pMC.endT  = 0.1;
pMC.nMC   = 1e3;
pMC.nORel = 2;

aStat = sq.anneal(pMC);

% energy after MC
E0 = sq.energy;

% further optimisation using steepest descendent method
sq.optmagsteep;
E1 = sq.energy;


%% plot magnetic structure

sq.plot('labelatom',false,'sspin',4)

%% spin waves

H = 1;
K = 1;
nQ = 200;

spec = sq.spinwave({[0 0 0] [H*N1 K*N2 0] nQ},'hermit',false);
spec = sw_egrid(spec);

figure
sw_plotspec(spec)


%% no disorder

sq2 = sw_model('squareAF',[1 0.3]);
spec = sq2.spinwave({[0 0 0] [H K 0] nQ});
spec = sw_egrid(spec);

%% plot spec

figure
sw_plotspec(spec,'mode','disp','imag',true,'colorbar',false,'axLim',[0 3],'colormap',[0 0 0])









