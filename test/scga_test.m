%% define the pyrochlore lattice

pyro = spinw;
pyro.genlattice('lat_const',3*sqrt(8)*[1 1 1])

% basis vectors of the FCC lattice
a = [0 0 0;0 1/2 1/2;1/2 0 1/2;1/2 1/2 0];


subIdx = repmat(1:4,1,4);
% colors of the sublattices
color  = {'red' 'green' 'orange' 'yellow'};
for ii = 1:4
    for jj = 1:4
        pyro.addatom('r',a(ii,:)+a(jj,:)/2,'S',1,'color',color{ii})
        
    end
end

pyro.quickham(1)

plot(pyro)

%% Fourier transform of the interaction matrix

nMag = numel(pyro.matom.idx);
idx1 = repmat(subIdx',1,nMag);
idx2 = repmat(subIdx,nMag,1);

Q = [1/3 1/4 1/5];
FT = pyro.fourier(Q');
% sum up the FT per sublattice

FT = reshape(FT.ft,9,[]);

subs = [idx1(:) idx2(:)];

FT2 = zeros(4,4,9);

for ii = 1:9
    FT2(:,:,ii) = accumarray(subs,FT(ii,:));
end

FT2 = sum(FT2,3)/12+eye(4);

cxy = cos((Q(1)+Q(2))/4*2*pi);
cxz = cos((Q(1)+Q(3))/4*2*pi);
cyz = cos((Q(2)+Q(3))/4*2*pi);

cxym = cos((Q(1)-Q(2))/4*2*pi);
cxzm = cos((Q(1)-Q(3))/4*2*pi);
cyzm = cos((Q(2)-Q(3))/4*2*pi);

FT3 = [1 cyz cxz cxy;cyz 1 cxym cxzm;cxz cxym 1 cyzm;cxy cxzm cyzm 1];

FT2-FT3

%% high temperature solution

kbT = sw_converter(1,'K','meV');
beta = 1/kbT;

lambda = 3;

%% solution of self-consistent equation

nQ0 = 1e4;
D  = 3;
% FT
N  = round(nQ0^(1/D));
nQ = N^D;
BZ = sw_qgrid('mat',eye(3),'bin',repmat({linspace(0,1,N)},1,D));

FT = pyro.fourier(reshape(BZ,3,[]));
FT = FT.ft;
% include the spin value into the Fourier transform of the Js
% thus convert the model into interacting S=1 spins
FT = bsxfun(@times,FT,permute(bsxfun(@times,pyro.matom.S',pyro.matom.S),[3 4 1 2]));


%%
subIdx = 1:16;
%subIdx = repmat(1:4,1,4);
% reduce the lattice into sublattices
nMag = numel(pyro.matom.idx);
idx1 = repmat(subIdx',1,nMag);
idx2 = repmat(subIdx,nMag,1);
subs = [idx1(:) idx2(:)];

nSub = max(subIdx);

FT  = reshape(FT,3,3,[],nQ);
FT2 = zeros(3,3,nSub,nSub,nQ);

for ii = 1:nSub
    for jj = 1:nSub
        FT2(:,:,ii,jj,:) = permute(sum(FT(:,:,ismember(subs,[ii jj],'rows'),:),3),[1 2 3 5 4])*nSub/nMag;
        if ii == jj
            FT2(:,:,ii,ii,:) = FT2(:,:,ii,ii,:) + 1;
        end
    end
end

[E,D] = eigorth(squeeze(FT2(1,1,:,:,:)));

%% plot the fitted curve

beta = 200;
lv = linspace(1,4,101);
for ii = 1:numel(lv)
    lambda = lv(ii);
    J      = 1;
    Di(ii) = abs(sumn(1./(lambda+beta*J*D),[1 2])/nQ/4-1/3);
    Di2(ii) = sumn(1./(lambda+beta*J*D),[1 2])/nQ/4;
end
fminsearch(@(lambda)abs(sumn(1./(lambda+beta*J*D/4),[1 2])/nQ/4-1/3),3)
figure
plot(lv,Di2)
hold on
plot(lv,lv*0+1/3)

%% find the values of lambda

Jv = 10.^linspace(-2,2,51);
beta = 1;

for ii = 1:numel(Jv)
    J = Jv(ii);
    li(ii) = fminsearch(@(lambda)abs(sumn(1./(lambda+beta*J*D),[1 2])/nQ/4-1/3),3);
end

figure;
semilogx(1./Jv,li,'o-')
xlabel('T/J_1')
ylabel('\lambda')
axis([1e-2 1e2 1 3.5])

%% calculate the spin correlations

%subIdx = 1:16;
subIdx = repmat(1:4,1,4);

beta = 1;
J    = 200;

lambda = fminsearch(@(lambda)abs(sumn(1./(lambda+beta*J*D),[1 2])/size(D,2)/4-1/3),3);

Q = sw_qgrid('bin',{[0 0.02 4] 0 [-1 0.02 4]});
%Q = sw_qgrid('u',[1 1 0],'v',[0 0 1],'bin',{[0 0.02 3] [-1 0.02 4]});
nQ = numel(Q)/3;

FT = pyro.fourier(reshape(Q,3,[]));
FT = FT.ft;
% include the spin value into the Fourier transform of the Js
% thus convert the model into interacting S=1 spins
FT = bsxfun(@times,FT,permute(bsxfun(@times,pyro.matom.S',pyro.matom.S),[3 4 1 2]));

% reduce the lattice into sublattices
nMag = numel(pyro.matom.idx);
idx1 = repmat(subIdx',1,nMag);
idx2 = repmat(subIdx,nMag,1);
subs = [idx1(:) idx2(:)];

nSub = max(subIdx);

FT  = reshape(FT,3,3,[],nQ);
FT2 = zeros(3,3,nSub,nSub,nQ);

for ii = 1:nSub
    for jj = 1:nSub
        FT2(:,:,ii,jj,:) = permute(sum(FT(:,:,ismember(subs,[ii jj],'rows'),:),3),[1 2 3 5 4])*nSub/nMag;
        if ii == jj
            FT2(:,:,ii,ii,:) = FT2(:,:,ii,ii,:) + 1;
        end
    end
end
V = FT2;

V = bsxfun(@plus,permute(lambda*eye(nSub),[3 4 1 2]),beta*J*V);

Sab = zeros(1,nQ);
for ii = 1:nQ
    Sab(ii) = sumn(inv(squeeze(V(1,1,:,:,ii))),[1 2]);
end

dQ = num2cell(size(Q));
% Correlations per site
Sab = reshape(Sab,dQ{2:end})/nSub;

%% plot correlations
clf
hSurf = surf(squeeze(Q(1,:,:)),squeeze(Q(3,:,:)),squeeze(Sab)*4);
hold on
contour3(squeeze(Q(1,:,:)),squeeze(Q(3,:,:)),squeeze(Sab)*4,0.5:0.5:2.5,'color','k');
hSurf.EdgeAlpha = 0;
view(2)
%axis([0 3 -1 4])
box on
colormap(cm_viridis)
caxis([0 2.7])
colorbar

%% plot cuts

figure
nCut = 4;
idxv = round(linspace(1,251,nCut));
for ii = 1:nCut
    plot(squeeze(Q(1,:,1,1)),Sab(:,1,idxv(ii))*4)
    hold on
end

%% test the new SCGA code

pyro.setunit('mode','1')

%Q = sw_qgrid('bin',{[0 0.02 4] 0 [-1 0.02 4]});
Q = sw_qgrid('u',[1 1 0],'v',[0 0 1],'bin',{[0 0.02 3] [-1 0.02 4]});

betaJ = 100;
%betaJ = 10.^linspace(-5,5,31);
subLat = repmat(1:4,1,4);
%subLat = 1:16;
tic
spec2 = pyro.scga2(Q,'T',1/2./betaJ,'plot',true,'nInt',1e4,'subLat',subLat,'lambda',[],'isomode','off');
toc
ylim([0 5])

%% plot lambda values

figure;
semilogx(spec.T*2,spec.lambda,'o-')
%xlim([1e-2 1e2])

%% plot correlations

spec = spec2;
Q   = spec.hkl;
Sab = squeeze(spec.Sab(1,1,:,:)*4);

figure
hSurf = surf(squeeze(Q(1,:,:)),squeeze(Q(3,:,:)),Sab);
hold on
contour3(squeeze(Q(1,:,:)),squeeze(Q(3,:,:)),Sab,0.5:0.5:2.5,'color','k');
hSurf.EdgeAlpha = 0;
view(2)
%axis([0 3 -1 4])
box on
colormap(cm_viridis)
caxis([0 2.7])
colorbar

%% test sublattice problem

%subLat = 1:16;
subLat = [];
%subLat = repmat(1:4,1,4);

beta = 200;

nMag = numel(pyro.matom.idx);
idx1 = repmat(subLat',1,nMag);
idx2 = repmat(subLat,nMag,1);
subs = [idx1(:) idx2(:)];

if ~isempty(subLat)
    nSub = max(subLat);
else
    nSub = nMag;
end

BZ   = sw_qgrid('bin',{[0 0.1 2] [0 0.1 2] 0.3});
sQBZ = num2cell(size(BZ));
nQBZ = numel(BZ)/3;
chi0 = pyro.fourier(reshape(BZ,3,[]),'extend',false);
% include the spin value into the Fourier transform of the Js
% thus convert the model into interacting S=1 spins
%FT = bsxfun(@times,chi0.ft,permute(bsxfun(@times,S',S),[3 4 1 2]));
FT = chi0.ft;

% reduce the lattice into sublattices
if ~isempty(subLat)
    FT  = reshape(FT,3,3,[],nQBZ);
    FT2 = zeros(3,3,nSub,nSub,nQBZ);
    
    for ii = 1:nSub
        for jj = 1:nSub
            FT2(:,:,ii,jj,:) = permute(sum(FT(:,:,ismember(subs,[ii jj],'rows'),:),3),[1 2 3 5 4])*nSub/nMag;
            if ii == jj
                FT2(:,:,ii,ii,:) = FT2(:,:,ii,ii,:) + 1;
            end
        end
    end
    FT = FT2;
else
    for ii = 1:nSub
        FT(:,:,ii,ii,:) = FT(:,:,ii,ii,:) + 1;
    end
end
% find the eigenvalues over the BZ
FT = squeeze(FT(1,1,:,:,:));
omega = zeros(nSub,nQBZ);
for ii = 1:nQBZ
    omega(:,ii) = eig(FT(:,:,ii));
end

% find the optimum value of lambda
lambda = fminsearch(@(lambda)abs(sumn(1./(lambda+beta*omega),[1 2])/nQBZ/nSub-1/3),3)

% plot omega
omega = reshape(omega',sQBZ{2:end},nSub);

%% simple eig

BZ   = sw_qgrid('bin',{[1 0.1 3] [1 0.1 3] 0.25});
sQBZ = num2cell(size(BZ));
nQBZ = numel(BZ)/3;
chi0 = pyro.fourier(reshape(BZ,3,[]),'extend',false);
% include the spin value into the Fourier transform of the Js
% thus convert the model into interacting S=1 spins
%FT = bsxfun(@times,chi0.ft,permute(bsxfun(@times,S',S),[3 4 1 2]));
FT = chi0.ft;
FT = squeeze(FT(1,1,:,:,:));

omega = zeros(nSub,nQBZ);
for ii = 1:nQBZ
    omega(:,ii) = eig(FT(:,:,ii));
end
% plot omega
omega = reshape(omega',sQBZ{2:end},nSub);

%%
clf
for ii = 1:nSub
    surf(squeeze(BZ(1,:,:)),squeeze(BZ(2,:,:)),omega(:,:,ii))
    hold on
end
%axis([0 2 0 2 -4 4])
ax=axis;
axis([ax(1:4) -4 4])
view(az,el)

%% test the new SCGA code

ybti = spinw('~/Documents/structures/Yb2Ti2O7/Yb2Ti2O7_cryst.cif');

ybti.setunit('mode','1')
ybti.quickham(1);
%Q = sw_qgrid('bin',{[0 0.02 4] 0 [-1 0.02 4]});
Q = sw_qgrid('u',[1 1 0],'v',[0 0 1],'bin',{[0 0.02 3] [-1 0.02 4]});

betaJ = 100;
%betaJ = 10.^linspace(-2,2,31);
%subLat = repmat(1:4,1,4);
subLat = 1:16;
spec = ybti.scga(Q,'T',1/2./betaJ,'plot',true,'nInt',1e4,'subLat',subLat,'lambda',[]);
ylim([0 5])
