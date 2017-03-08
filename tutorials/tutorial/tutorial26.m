% Spiral dispersion and correlation function on Bravais lattice
% Theoretical and experimental calculation

%% chain

J1 = 1;
J2 = 3;

ch = spinw;
ch.genlattice('lat_const',[3 4 4])
ch.addatom('r',[0 0 0],'S',1)
ch.addmatrix('label','J1','value',1)
ch.addmatrix('label','J2-','value',3,'color','blue')

ch.gencoupling
ch.addcoupling('mat','J1','bond',1)
ch.addcoupling('mat','J2-','bond',5)

plot(ch,'range',[3 1 1])

%%

ch.optmagstr('func',@gm_planar,'xmin',[0, 0 0 0, 0 0],'xmax',[0 1/2 0 0, 0 0])

%%

spec = ch.spinwave({[0 0 0] [1 0 0] 500});
spec = sw_egrid(spec);

figure
subplot(2,1,1)
sw_plotspec(spec,'mode','disp','linestyle','-');
subplot(2,1,2)
sw_plotspec(spec,'mode','int','linestyle','-');

%

Q = ch.mag_str.k(1);

J = @(k)2*J1*cos(2*pi*k)+2*J2*cos(4*pi*k);
A = @(k)J(k)+J(k+Q)/2+J(k-Q)/2-2*J(Q);
B = @(k)J(k)-J(k+Q)/2-J(k-Q)/2;

k = linspace(0,1,500);

w1 =  @(k)sqrt(A(k).^2-B(k).^2);
w2 =  @(k)-sqrt(A(k).^2-B(k).^2);
%w1 = A(k)+B(k);
%w2 = A(k)-B(k);
Sxx   = 2*(A(k)-B(k))./w1(k);
Sxixi = @(k)2*(A(k)+B(k))./w1(k);
Syy   = 1/4*(Sxixi(k-Q)+Sxixi(k+Q));
Szz   = Syy;

%figure;
subplot(2,1,1)
hold on
plot(k,(w1(k))/2,'ro')
hold on
plot(k,(w1(k-Q))/2,'go')
plot(k,(w1(k+Q))/2,'bo')
axis([0 1 0 10])
%
subplot(2,1,2)
hold on
plot(k,Sxx/4+0.01,'r-')
hold on
plot(k,Sxixi(k-Q)/4/4+0.01,'g-')
plot(k,Sxixi(k+Q)/4/4+0.01,'b-')

axis([0 1 0 1])
