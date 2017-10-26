function spec = partest01(N)

if nargin == 0
    N = 1e3;
end

hPool = gcp('nocreate');
if isempty(hPool)
    nWorker = 1;
else
    nWorker = hPool.NumWorkers;
end

N = round(N/nWorker)*nWorker;

% setup
addpath(genpath('/home/lvi05884/spinw'))
swpref.setpref('usemex',false,'tid',0,'fid',0);

yig = yig_create;

Q = rand(3,N);

fprintf('spinwave N=%d, nomex\n',N)
tic
spec1 = yig.spinwave(Q(:,1));
toc

fprintf('spinwavefast N=%d, nomex\n',N)
tic
spec2 = yig.spinwavefast(Q);
toc

fprintf('spinwavefast_duc N=%d, nomex\n',N)
tic
spec3 = spinwavefast_duc(yig,Q);
toc

fprintf('spinwave N=%d, nomex, spmd\n',N)
tic

Qc = Composite();
Ni = N/nWorker;
for ii = 1:nWorker
    Qc{ii} = Q(:,(1:Ni)+(ii-1)*Ni);
end

spmd
    swpref.setpref('fid',0,'tid',0);
    specSpmd = yig.spinwave(Qc);
    Sab = specSpmd.Sab;
end

spec4.Sab = cat(4,Sab{:});
toc

spec = {spec1 spec2 spec3};

end

function yig = yig_create

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

end