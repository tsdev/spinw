%% kagome lattice

kag = spinw;
kag.genlattice('lat_const',[6 6 4],'angled',[90 90 120])
kag.addatom('r',[1/2 0   0],'S',1)
kag.addatom('r',[1/2 1/2 0],'S',1)
kag.addatom('r',[0   1/2 0],'S',1)

%kag.newcell('bvect',{[2 0 0] [0 1 0] [0 0 1]})
kag.quickham(1)
%plot(kag)

%% eig

Q = sw_qgrid('bin',{[0 0.05 2] [0 0.05 2]});

chi = kag.fourier(reshape(Q,3,[]));

ft = squeeze(chi.ft(1,1,:,:,:));

clear om
for ii = 1:size(ft,3)
    om(:,ii) = eig(ft(:,:,ii));
end

nMode = size(om,1);
sQ = num2cell(size(Q));
om = reshape(om',sQ{2:end},nMode);

clf
for ii = 1:nMode
    surf(squeeze(Q(1,:,:)),squeeze(Q(2,:,:)),om(:,:,ii));
    hold on
end
%% scga method determine lambda

kag.setunit('mode','1')
T = 10.^linspace(-2,2,41);
Q = sw_qgrid('bin',{[0 0.05 2] [0 0.05 2]});
spec = kag.scga(Q,'T',T,'plot',true,'nInt',1e4);

figure
semilogx(spec.T,spec.lambda,'o-')

%% scga method diffuse scattering

kag.setunit('mode','1')
T = 1;
Q = sw_qgrid('bin',{[0 0.05 4] [0 0.05 4]});
spec = kag.scga2(Q,'T',T,'plot',true,'nInt',1e4);


%%
figure
hSurf = surf(squeeze(Q(1,:,:,:)),squeeze(Q(2,:,:,:)),spec.Sab);
hSurf.EdgeAlpha = 0;
view(2)
hold on
contour3(squeeze(Q(1,:,:,:)),squeeze(Q(2,:,:,:)),spec.Sab,0:0.05:2,'color','k');




