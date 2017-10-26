%% yig setup

yig = spinw('https://goo.gl/kQO0FJ');
% spin quantum number of Fe3+ ions, determined automatically by SpinW
S0 = max(yig.unit_cell.S);
% normalize spins to S=1 as it is in the paper
yig.unit_cell.S = yig.unit_cell.S/S0;
% new basis vectors in rows
pBV = [1/2 1/2 -1/2;-1/2 1/2 1/2;1/2 -1/2 1/2];
% lattice constant of YIG
% Generate the bonds using centered cell
% An interesting symmetry property of YIG in the "I a -3 d" space group is
% that the bonds type 3 and type 4 have the exact same length, however they
% are not related by symmetry. This can be easily seen by checking the
% center psition of the bonds:
yig.gencoupling('maxDistance',6);
% Create spin Hamiltonian
% change from BCC to primitive cubic cell
T = yig.newcell('bvect',{pBV(1,:) pBV(2,:) pBV(3,:)});
% exchange values from the paper
Jad = sw_converter(9.60e-21,'J','THz','photon');
Jdd = sw_converter(3.24e-21,'J','THz','photon');
Jaa = sw_converter(0.92e-21,'J','THz','photon');

% scale the interactions from classical moment size to quantum model
Scl = sqrt(S0*(S0+1));
yig.quickham([Jad Jdd Jaa]/Scl)

% add external field and convert from the standard SpinW unit (meV) to THz
yig.field([0 0 0.01]*sw_converter(1,'meV','THz','photon'))
yig.optmagsteep('nRun',1e2)
yig.genmagstr('mode','rotate','n',[0 0 1])

if sum(yig.mag_str.F(3,:),2)<0
    yig.mag_str.F = -yig.mag_str.F;
end

yig.mag_str.F = real(yig.mag_str.F);
yig.mag_str.F(1:2,:) = 0;
yig.mag_str.F = bsxfun(@rdivide,yig.mag_str.F,abs(yig.mag_str.F(3,:)));

%%

Q0  = T*[1 2 3]';
Q_N = T*[ 1/2  1/2    0]'+Q0;
Q_G = T*[   0    0    0]'+Q0;
Q_H = T*[   0    0    1]'+Q0;
N = 1e3;
Q = {Q_N Q_G Q_H N};

swpref.setpref('usemex',true);

spec1 = yig.spinwave(Q);

