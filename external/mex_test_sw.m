% Test systems are:
%   ndlt811:        Laptop - Dell Precision M3800; 4xi7-4712HQ; 16GB, Windows 7
%   eryenyo:        Laptop - Gigabyte P34; 4xi7-4700HQ; 16GB, Ubuntu 14.10
%   ndl01wkc26243:  Workstation - Dell Precision T7500; 12xXeon-X5650; 50GB, CentOS 7

%setenv('OMP_NUM_THREADS','24'); clear mex;
%addpath(genpath('~/spinw'))
%mex('-g','-v','-largeArrayDims','eig_omp.cpp','-lmwlapack')
%mex('-g','-v','-largeArrayDims','chol_omp.cpp','-lmwlapack','-lmwblas')
%clear mex;
addpath(genpath('~/tmp/spinw'))

% Antiferromagnetic square lattice (tutorial 4, Cu+1 S=1) - small system [3 spins, n=6]
AFsq = sw; 
AFsq.genlattice('lat_const',[3 3 10],'angled',[90 90 90],'sym',0)
AFsq.addatom('r',[0 0 0],'S', 1,'label','Cu1','color','b');
AFsq.gencoupling('maxDistance',9)
AFsq.addmatrix('label','J1','value',1,'color','red');      AFsq.addcoupling('J1',1)
AFsq.addmatrix('label','J2','value',-0.1,'color','green'); AFsq.addcoupling('J2',2)
AFsq.genmagstr('mode','helical','k',[1/2 1/2 0],'n',[0 0 1], 'S',[1; 0; 0],'nExt',[1 1 1]);  
%plot(AFsq,'range',[2 2 0.5],'zoom',-1)
%AFsq.addmatrix('mat',rand(3),'label',{'g0'}); AFsq.addg('g0');

itst=1;
% Runs test
hkl = {[1/4 3/4 0] [1/2 1/2 0] [1/2 0 0] [3/4 1/4 0] [1 0 0] [3/2 0 0] 10000};
nm = 5;   % Ensure same number of slices for all tests
tic; linespec_herm_mex = AFsq.spinwave(hkl,'hermit',true,'useMex',true,'optmem',nm); t1=toc; linespec_herm_mex = sw_egrid(sw_neutron(linespec_herm_mex));
tic; linespec_herm_nomex = AFsq.spinwave(hkl,'hermit',true,'useMex',false,'optmem',nm); t2=toc; linespec_herm_nomex = sw_egrid(sw_neutron(linespec_herm_nomex));
tic; linespec_nonherm_mex = AFsq.spinwave(hkl,'hermit',false,'useMex',true,'optmem',nm); t3=toc; linespec_nonherm_mex = sw_egrid(sw_neutron(linespec_nonherm_mex));
tic; linespec_nonherm_nomex = AFsq.spinwave(hkl,'hermit',false,'useMex',false,'optmem',nm); t4=toc; linespec_nonherm_nomex = sw_egrid(sw_neutron(linespec_nonherm_nomex));
tt{itst}=[t1 t2 t3 t4]; 
ch1(itst,:)=[sum(abs(linespec_herm_mex.omega(:)-linespec_herm_nomex.omega(:))) sum(abs(linespec_nonherm_mex.omega(:))-abs(linespec_nonherm_nomex.omega(:)))];
ch2(itst,:)=[sum(abs(linespec_herm_mex.swConv(:)-linespec_herm_nomex.swConv(:))) sum(abs(linespec_nonherm_mex.swConv(:)-linespec_nonherm_nomex.swConv(:)))];

disp('---------------------------------------------------------------------------------------');
disp(sprintf('             %16s  %16s  %16s  %16s','Hermitian Mex','Hermitian NoMex','NonHermitian Mex','NonHermitian NoMex')); 
disp(sprintf('Run Time(s)  % 16.6f  % 16.6f  % 16.6f    % 16.6f',tt{itst}));
disp(sprintf('Check omega  %13s     % 16.6g  %13s       % 16.6g',' ',ch1(itst,1),' ',ch1(itst,2)));
disp(sprintf('Check intensities %13s% 16.6g  %13s       % 16.6g',' ',ch2(itst,1),' ',ch2(itst,2)));
disp('---------------------------------------------------------------------------------------');

% Times
%   ndlt811:         43.7453   72.0950  406.1366  611.1799
%   eryenyo:         24.1470   85.5668  382.7930  563.6103
%   ndl01wkc26243:   36.2106  107.8983  527.3301  796.2985

% KCu3As2O7(OD)3 kagome (Nilsen PRB 89 140412) - Tutorial 18, medium sized [18 spins, n=36]
J   = -2; Jp  = -1; Jab = 0.75; Ja  = -J/.66 - Jab; Jip = 0.01;
hK = sw;
hK.genlattice('lat_const',[10.2 5.94 7.81],'angled',[90 117.7 90],'sym','C 2/m');
hK.addatom('r',[0   0   0],'S',1/2,'label','MCu2','color','b');
hK.addatom('r',[1/4 1/4 0],'S',1/2,'label','MCu2','color','k');
hK.gencoupling;
hK.addmatrix('label','J',  'color','r',   'value',J);   hK.addcoupling('J',1);
hK.addmatrix('label','J''','color','g',   'value',Jp);  hK.addcoupling('J''',2);
hK.addmatrix('label','Ja', 'color','b',   'value',Ja);  hK.addcoupling('Ja',3);
hK.addmatrix('label','Jab','color','cyan','value',Jab); hK.addcoupling('Jab',5);
hK.addmatrix('label','Jip','color','gray','value',Jip); hK.addcoupling('Jip',10);
hK.genmagstr('mode','helical','n',[0 0 1],'S',[1 0 0]','k',[0.77 0 0.115],'next',[1 1 1]);
optpar.func = @gm_planar;
optpar.nRun = 5;
optpar.xmin = [    zeros(1,6), 0.5 0 0.0, 0 0];
optpar.xmax = [2*pi*ones(1,6), 1.0 0 0.5, 0 0];
magoptOut = hK.optmagstr(optpar);
kOpt = hK.mag_str.k;
hK.genmagstr('mode','helical','n',[0 0 1],'S',[1 0 0]','k',kOpt,'next',[1 1 1]);
%plot(hK,'range',[2 2 0.3],'sSpin',2)

itst=2;
% Runs test
hkl = {[0 0 0] [1 0 0] 50000};
nm = 10;  % Ensure same number of slices for all tests

tic; linespec_herm_mex = hK.spinwave(hkl,'hermit',true,'useMex',true,'optmem',nm); t1=toc; linespec_herm_mex = sw_egrid(sw_neutron(linespec_herm_mex));
tic; linespec_herm_nomex = hK.spinwave(hkl,'hermit',true,'useMex',false,'optmem',nm); t2=toc; linespec_herm_nomex = sw_egrid(sw_neutron(linespec_herm_nomex));
tic; linespec_nonherm_mex = hK.spinwave(hkl,'hermit',false,'useMex',true,'optmem',nm); t3=toc; linespec_nonherm_mex = sw_egrid(sw_neutron(linespec_nonherm_mex));
tic; linespec_nonherm_nomex = hK.spinwave(hkl,'hermit',false,'useMex',false,'optmem',nm); t4=toc; linespec_nonherm_nomex = sw_egrid(sw_neutron(linespec_nonherm_nomex));
tt{itst}=[t1 t2 t3 t4]; 
ch1(itst,:)=[sum(abs(linespec_herm_mex.omega(:)-linespec_herm_nomex.omega(:))) sum(abs(linespec_nonherm_mex.omega(:))-abs(linespec_nonherm_nomex.omega(:)))];
ch2(itst,:)=[sum(abs(linespec_herm_mex.swConv(:)-linespec_herm_nomex.swConv(:))) sum(abs(linespec_nonherm_mex.swConv(:)-linespec_nonherm_nomex.swConv(:)))];

disp('---------------------------------------------------------------------------------------');
disp(sprintf('             %16s  %16s  %16s  %16s','Hermitian Mex','Hermitian NoMex','NonHermitian Mex','NonHermitian NoMex')); 
disp(sprintf('Run Time(s)  % 16.6f  % 16.6f  % 16.6f    % 16.6f',tt{itst}));
disp(sprintf('Check omega  %13s     % 16.6g  %13s       % 16.6g',' ',ch1(itst,1),' ',ch1(itst,2)));
disp(sprintf('Check intensities %13s% 16.6g  %13s       % 16.6g',' ',ch2(itst,1),' ',ch2(itst,2)));
disp('---------------------------------------------------------------------------------------');

%{
[~,nSuperlat] = rat(hK.mag_str.k,0.01);
hK.genmagstr('mode','helical','next',nSuperlat)
hK.mag_str.k = [0 0 0];
out = optmagsteep(hK,'nRun',200); hK = out.obj;

hkl = {[0 0 0] [1 0 0] 40};
nm = 10;  % Ensure same number of slices for all tests
tic; linespec_herm_mex = hK.spinwave(hkl,'hermit',true,'useMex',true,'optmem',nm); t5=toc; 
tic; linespec_herm_nomex = hK.spinwave(hkl,'hermit',true,'useMex',false,'optmem',nm); t6=toc; 
tic; linespec_nonherm_mex = hK.spinwave(hkl,'hermit',false,'useMex',true,'optmem',nm); t7=toc; 
tic; linespec_nonherm_nomex = hK.spinwave(hkl,'hermit',false,'useMex',false,'optmem',nm); t8=toc; 
disp(sprintf('Supercell    % 16.6f  % 16.6f  % 16.6f    % 16.6f',t5,t6,t7,t8)); 
%}

% Times
%   ndlt811:         26.4661   40.5760   60.6112   63.3341
%   eryenyo:          7.6710   36.0402   56.3213   73.5568
%   ndl01wkc26243:   17.2135   42.8147   74.9062  106.8541

% Bi4Fe5O13F - large(ish) system [80 spins, n=160]
ff=11.6*2.5; Jc1 = 34/ff; Jc2 = 20/ff; Jab1 = 45/ff; Jab2 = 74/ff; Jd = 191/ff;
bfof = sw();
bfof.genlattice('lat_const',[8.29950 8.29950 18.05730],'angled',[90 90 90],'sym','P 42/m b c');
bfof.addatom('r',[0.5  0.   0.0800],'S',2.5,'color',[0 0 255],'label','Fe1');
bfof.addatom('r',[0.8515 0.8388 0], 'S',2.5,'color',[255 0 0],'label','Fe2');
bfof.addatom('r',[0.5  0.   0.25  ],'S',2.5,'color',[0 0 128],'label','Fe1_3');
bfof.gencoupling('maxBond',99,'maxDistance',10);
bfof.addmatrix('value',Jc1,'label','Jc1','color','r');         bfof.addcoupling('Jc1',1);
bfof.addmatrix('value',Jc2,'label','Jc2','color',[128 0 0]);   bfof.addcoupling('Jc2',2);
bfof.addmatrix('value',Jab1,'label','Jab1','color','b');       bfof.addcoupling('Jab1',3);
bfof.addmatrix('value',Jab2,'label','Jab2','color',[0 255 0]); bfof.addcoupling('Jab2',4);
bfof.addmatrix('value',Jd,'label','Jd','color','k');           bfof.addcoupling('Jd',5);
bfof.addmatrix('value',diag([0 0 0.2]),'label','D');           bfof.addaniso('D');
S2a = -[4.05 -0.35 0]; S2b = -[-0.35 4.05 0]; S1a = [2.18 -2.53 0]; S1b = [2.53 2.18 0];
S = [S1b; S1a; S1a; S1b; S1b; S1a; S1a; S1b; S2a; -S2a; -S2b; S2b; S2b; -S2b; S2a; -S2a; -S1b; -S1a; -S1b; -S1a];
Sv=[S; -S; -S; S];
bfof.genmagstr('mode','direct','S',Sv','nExt',[2 2 1]);
%bfof.plot('range',[0 2; 0 2; 0.0 1.0]); view([45 90]); set(gcf,'Tag','');
out = optmagsteep(bfof,'nRun',200);
%plot(out.obj,'range',[0 2; 0 2; 0.0 1.0]); view([45 90]);
bfof = out.obj;

itst=3;
% Runs test
nm = 6;   % For the workstation with 50GB.
nm = 60;  % For laptops with 16GB memory, to have the same number of slices
hkl={[0 0 0] [1 1 0] [1 1 1] [0 0 1] 2000};
tic; linespec_herm_mex = bfof.spinwave(hkl,'hermit',true,'useMex',true,'optmem',nm); t1=toc; linespec_herm_mex = sw_egrid(sw_neutron(linespec_herm_mex)); 
tic; linespec_herm_nomex = bfof.spinwave(hkl,'hermit',true,'useMex',false,'optmem',nm); t2=toc; linespec_herm_nomex = sw_egrid(sw_neutron(linespec_herm_nomex)); 
tic; linespec_nonherm_mex = bfof.spinwave(hkl,'hermit',false,'useMex',true,'optmem',nm); t3=toc; linespec_nonherm_mex = sw_egrid(sw_neutron(linespec_nonherm_mex)); 
tic; linespec_nonherm_nomex = bfof.spinwave(hkl,'hermit',false,'useMex',false,'optmem',nm); t4=toc; linespec_nonherm_nomex = sw_egrid(sw_neutron(linespec_nonherm_nomex)); 
tt{itst}=[t1 t2 t3 t4]; 
ch1(itst,:)=[sum(abs(linespec_herm_mex.omega(:)-linespec_herm_nomex.omega(:))) sum(abs(linespec_nonherm_mex.omega(:))-abs(linespec_nonherm_nomex.omega(:)))];
ch2(itst,:)=[sum(abs(linespec_herm_mex.swConv(:)-linespec_herm_nomex.swConv(:))) sum(abs(linespec_nonherm_mex.swConv(:)-linespec_nonherm_nomex.swConv(:)))];

for itst=1:3
disp('---------------------------------------------------------------------------------------');
disp(sprintf('             %16s  %16s  %16s  %16s','Hermitian Mex','Hermitian NoMex','NonHermitian Mex','NonHermitian NoMex')); 
disp(sprintf('Run Time(s)  % 16.6f  % 16.6f  % 16.6f    % 16.6f',tt{itst}));
disp(sprintf('Check omega  %13s     % 16.6g  %13s       % 16.6g',' ',ch1(itst,1),' ',ch1(itst,2)));
disp(sprintf('Check intensities %13s% 16.6g  %13s       % 16.6g',' ',ch2(itst,1),' ',ch2(itst,2)));
disp('---------------------------------------------------------------------------------------');
end

% Times
%   ndlt811:        123.6812  139.6476  197.1951  363.2818
%   eryenyo:        127.4535  135.8621  220.9112  389.4686
%   ndl01wkc26243:   67.8451  117.7316  146.2180  474.1894
