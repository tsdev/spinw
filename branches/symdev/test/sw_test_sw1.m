function [Res, errMsg] = sw_test_sw1(tol)
% check spin wave calculation (spinwave routine)

Res = 0;
errMsg = {};

try
    
    %% define crystal structure of Na2IrO3
    
    % data from: Choi, S. K., Coldea, et al. (2012). PRL, 108(12), 127204
    % space group C2/m, crystallographic parameters at 300 K
    
    nairo = sw;
    nairo.genlattice('lat_const',[5.427 9.395 5.614],'angle',[90 109.037 90]*pi/180,'sym','C 2/m');
    
    % add magnetic Ir
    nairo.addatom('r',[1/2; 0.167; 0],'S',1/2,'label','Ir','color',[140; 28; 22]);
    nairo.addatom('r',[0 1/2 1/2;0 0 0.340; 0 1/2 1/2],'S',[0 0 0],'label',{'Na1' 'Na2' 'Na3'},'color',ones(3)*233);
    nairo.addatom('r',[0.748 0.711; 0.178 0; 0.789 0.204],'S',[0 0],'label',{'O1', 'O2'},'color',[45 45;118 118;125 125]);
    
    hFig = plot(nairo);
    close(hFig);
    
    % define prototype magnetic Hamiltonian
    
    % regenerate crystal with P1 symmetry since the Hamiltonian is incompatible
    % with symmetry
    nairo.newcell({[1 0 0] [0 1 0] [0 0 1]})
    
    % rotation needs to be added, see sw_rot
    nairo.addmatrix('color',[255; 0; 0],'label','JKxx');
    nairo.addmatrix('color',[0; 255; 0],'label','JKyy');
    nairo.addmatrix('color',[0; 0; 255],'label','JKzz');
    
    % Heisenberg terms
    nairo.addmatrix('color',[128; 128; 128],'label','J1');
    nairo.addmatrix('color',[200; 128; 128],'label','J2');
    nairo.addmatrix('color',[128; 200; 128],'label','J3');
    
    % generate couplings up to 8 Angstrom
    nairo.gencoupling('maxdistance',8);
    
    % add J1, J2 and J3 and JK couplings
    nairo.addcoupling('J1',[1 2]);
    nairo.addcoupling('J2',[3 4]);
    nairo.addcoupling('J3',[7 8]);
    
    % anisotropic Kitaev couplings, incompatible with crystal symmetry
    nairo.addcoupling('JKxx',1,[1 4]);
    nairo.addcoupling('JKyy',1,[2 3]);
    nairo.addcoupling('JKzz',2);
    
    hFig = plot(nairo,'range',[-0.1 3.1;-0.1 3.1;-0.1 0.1]);
    close(hFig);
    
    % define Q scans
    nQ = 200;
    nE = 800;
    dQ = 0.764e-6;
    Qp{1} = [ -1;   0; 0]+dQ;
    Qp{2} = [  0;   0; 0]+dQ;
    Qp{3} = [  0;   1; 0]+dQ;
    Qp{4} = [  1;   1; 0]+dQ;
    Qp{5} = [1/2; 1/2; 0]+dQ;
    Qp{6} = [  0;   0; 0]+dQ;
    
    % case d stripy order
    % energy per spin is -0.2912
    J1 =  1; J2 =  0; J3 =  0; JK =  1.33;
    Na2IrO3fun(nairo,[J1 J2 J3 JK]);
    
    opt_par.xmin = [0 0 0 0,        0 0 0, 0 0];
    opt_par.xmax = [[1 1 1 1]*2*pi, 0 0 0, pi 0];
    opt_par.func = @gm_planar;
    opt_par.nRun = 10;
    
    nairo.optmagstr(opt_par);
    
    % calculate spin waves with zig-zag order
    E  = linspace(0,4,nE);
    
    specD = nairo.spinwave([Qp {nQ}]);
    specD = sw_neutron(specD,'pol',false);
    specD = sw_conv(specD,'convmode','Sxx','evect',E);
    
    hFig = figure;
    sw_plotspec(specD,'mode',3,'ahandle',gca,'imag',false,'convE',0.05,'axLim',0.5);
    sw_plotspec(specD,'mode',1,'ahandle',gca,'imag',false,'dashed',true,'colorbar',false,'legend',false,'title',false);
    close(hFig);
    
    %% case e-f stripy order
    % energy per spin is -0.375
    
    J1 =  1; J2 =  0; J3 =  0; JK =  2;
    Na2IrO3fun(nairo,[J1 J2 J3 JK]);
    
    opt_par.xmin = [0 0 0 0,        0 0 0, 0 0];
    opt_par.xmax = [[1 1 1 1]*2*pi, 0 0 0, pi 0];
    opt_par.func = @gm_planar;
    opt_par.nRun = 10;
    
    nairo.optmagstr(opt_par);
    
    % calculate spin waves with zig-zag order
    E  = linspace(0,4,nE);
    
    specEF = nairo.spinwave([Qp {nQ}]);
    specEF = sw_neutron(specEF,'pol',false);
    specEF = sw_conv(specEF,'convmode','Syy','evect',E);
    
    hFig = figure;
    sw_plotspec(specEF,'mode',3,'ahandle',gca,'imag',false,'convE',0.05,'axLim',0.5);
    sw_plotspec(specEF,'mode',1,'ahandle',gca,'imag',false,'dashed',true,'colorbar',false,'legend',false,'title',false);
    close(hFig);
    
    
    %% case g stripy order
    % optimise magnetic structure assuming planar
    % --> finds tripy order
    
    J1 =  1; J2 =  0.26; J3 = -0.2; JK = 0;
    Na2IrO3fun(nairo,[J1 J2 J3 JK]);
    
    opt_par.xmin = [0 0 0 0,        0 0 0, 0 0];
    opt_par.xmax = [[0 1 1 1]*2*pi, 0 0 0, 0 0];
    opt_par.func = @gm_planar;
    opt_par.nRun = 10;
    
    nairo.optmagstr(opt_par);
    
    % calculate spin waves with stripy order
    E  = linspace(0,2,nE);
    
    specG = nairo.spinwave([Qp {nQ}]);
    specG = sw_neutron(specG,'pol',false);
    specG = sw_conv(specG,'convmode','Szz','evect',E);
    
    hFig = figure;
    sw_plotspec(specG,'mode',3,'ahandle',gca,'imag',false,'convE',0.05,'axLim',0.5);
    sw_plotspec(specG,'mode',1,'ahandle',gca,'imag',false,'dashed',true,'colorbar',false,'legend',false,'title',false);
    close(hFig);
    
    %% case h zig-zag order
    % optimise magnetic structure assuming planar
    % --> finds tripy order
    
    J1 =  1; J2 =  0.78; J3 =  0.9; JK =  0;
    Na2IrO3fun(nairo,[J1 J2 J3 JK]);
    
    opt_par.xmin = [0 0 0 0,        0 0 0, 0 0];
    opt_par.xmax = [[0 1 1 1]*2*pi, 0 0 0, 0 0];
    opt_par.func = @gm_planar;
    opt_par.nRun = 10;
    
    nairo.optmagstr(opt_par);
    
    % calculate spin waves with zig-zag order
    E  = linspace(0,4,nE);
    
    specH = nairo.spinwave([Qp {nQ}]);
    specH = sw_neutron(specH,'pol',false);
    specH = sw_conv(specH,'convmode','Syy','evect',E);
    
    hFig = figure;
    subplot(2,1,1)
    sw_plotspec(specH,'mode',3,'ahandle',gca,'imag',false,'convE',0.05,'axLim',0.5);
    sw_plotspec(specH,'mode',1,'ahandle',gca,'imag',false,'dashed',true,'colorbar',false,'legend',false,'title',false,'colormap',[0 0 0]);
    
    % plot exact solution below
    subplot(2,1,2);
    specSimH = specH;
    specSimH.omega = repmat(omegaH([Qp {nQ}],[J1 J2 J3 JK]),4,1);
    sw_plotspec(specSimH,'mode',1,'ahandle',gca,'imag',false,'dashed',true,'colorbar',false,'colormap',[255 0 0]);
    title('Exact spin wave dispersion','fontsize',14);
    
    close(hFig);
    
    specH    = sw_omegasum(specH);
    specSimH = sw_omegasum(specSimH);
    
    % calculate difference between exact and numerical solutions
    omegaDiff = abs(specH.omega-specSimH.omega);
    omegaDiff(isnan(omegaDiff)) = 0;
    
    omegaDiff = max(max(omegaDiff));
    ratioH = omegaDiff/max(max(abs(real(specH.omega))));
    
    %% powder spectra
    J1 =  4.17; J2 =  0.78*J1; J3 =  0.9*J1; JK =  0;
    Na2IrO3fun(nairo,[J1 J2 J3 JK]);
    
    hklA = linspace(0.2,1.5,30);
    
    specHpow = nairo.powspec(hklA,'evect',linspace(0,8,100),'nRand',1e3);
    
    hFig = figure;
    sw_plotspec(specHpow)
    close(hFig);
    
    %% case ij zig-zag order
    % optimise magnetic structure assuming planar
    
    % DISPERSION DEPENDS ON THE SPIN DIRECTION
    % THE M POINT IS NOT UNIQUE
    
    J1 =  1; J2 =  0.23; J3 =  0.51; JK =  1.33;
    Na2IrO3fun(nairo,[J1 J2 J3 JK]);
    
    opt_par.xmin = [0 0 0 0,        0 0 0, 0 0];
    opt_par.xmax = [[0 1 1 1]*2*pi, 0 0 0, pi 0];
    opt_par.func = @gm_planar;
    opt_par.nRun = 10;
    
    nairo.optmagstr(opt_par);
    
    % calculate spin waves with zig-zag order
    E  = linspace(0,3,nE);
    
    specIJ = nairo.spinwave([Qp {nQ}]);
    specIJ = sw_neutron(specIJ,'pol',false);
    
    convmode = {'Sxx' 'Szz'};
    
    hFig = figure;
    
    for ii = 1:2
        subplot(3,1,ii);
        specIJ = sw_conv(specIJ,'convmode',convmode{ii},'evect',E);
        
        sw_plotspec(specIJ,'mode',3,'ahandle',gca,'imag',false,'convE',0.05,'axLim',0.5);
        sw_plotspec(specIJ,'mode',1,'ahandle',gca,'imag',false,'dashed',true,'colorbar',false,'legend',false,'title',false,'colormap',[0 0 0]);
        caxis([0 0.5])
    end
    
    % plot exact solution below
    % EXACT SOLUTION HAS WRONG EQUATIONS IN PUBLICATION
    subplot(3,1,3);
    specSimIJ = specIJ;
    omegaTemp = repmat(omegaIJ([Qp {nQ}],[J1 J2 J3 JK]),2,1);
    omegaTemp = bsxfun(@times,omegaTemp,[1 1 1 1 0 0 0 0]');
    specSimIJ.omega = omegaTemp;
    
    sw_plotspec(specSimIJ,'mode',1,'ahandle',gca,'imag',false,'dashed',true,'colorbar',false,'colormap',[255 0 0]);
    title('Exact spin wave dispersion','fontsize',14);
    
    close(hFig);
    
    % specIJ    = sw_omegasum(specIJ);
    % specSimIJ = sw_omegasum(specSimIJ);
    %
    % omegaDiff = abs(specIJ.omega-specSimIJ.omega);
    % omegaDiff(isnan(omegaDiff)) = 0;
    %
    % omegaDiff = max(max(omegaDiff));
    % ratioIJ = omegaDiff/max(max(abs(real(specIJ.omega))));
catch errMsg
    % code throws error
    Res = 1;
    return;
end

try
    if abs(ratioH) > tol
        error('sw_test_sw1:DataError','The calculated spin wave dispersion differs from the exact result!');
    end
catch errMsg
    Res = 2;
    return;
end

end


function Na2IrO3fun(obj, Jinp)
% takes the values as [J1 J2 J3 JK] and put them into the Hamiltonian

J1 = Jinp(1);
J2 = Jinp(2);
J3 = Jinp(3);
JK = Jinp(4);

JKxx=zeros(3); JKxx(1,1)=-JK;
JKyy=zeros(3); JKyy(2,2)=-JK;
JKzz=zeros(3); JKzz(3,3)=-JK;

obj.matrix.mat = cat(3,JKxx,JKyy,JKzz,eye(3)*J1,eye(3)*J2,eye(3)*J3);

end


function omega = omegaH(q,x)

q = sw_qscan(q);

J1 = x(1);
J2 = x(2);
J3 = x(3);

S = 1/2;

h = q(1,:);
k = q(2,:);

eta = exp(k*pi*1i/3);

A = S*(-J1 + 2*J2 + 3*J3 + 2*J2*cos(2*pi*h));
B = 2*S*J1*eta.*cos(pi*h);
C = 2*S*J2*(cos(pi*(h+k))+cos(pi*(h-k)));
D = S*(J1*eta.^2 + J3*(eta.^(-4)+2*eta.^2.*cos(2*pi*h)));

sum1 = A.^2 + B.*conj(B) - C.^2 - D.*conj(D);
sum2 = sqrt(4*abs(A.*B-C.*conj(D)).^2-abs(conj(B).*conj(D)-B.*D).^2);

omega = [sqrt(sum1 + sum2); sqrt(sum1 - sum2)];

end

function omega = omegaIJ(q,x)

q = sw_qscan(q);

J1 = x(1);
J2 = x(2);
J3 = x(3);
JK = x(4);

S = 1/2;

h = q(1,:);
k = q(2,:);

eta = exp(k*pi*1i/3);
xi  = exp(h*pi*1i);

A = S*(-J1 + 2*J2 + 3*J3 + 2*J2*cos(2*pi*h) + JK);
B = -1/2*S*JK*eta.^(-2);
C = S*(2*J1*cos(pi*h)-1/2*JK*xi).*eta;
D = S*(J1*eta.^(-2) + J3*(eta.^4+2*eta.^(-2).*cos(2*pi*h))-JK/2*eta.^(-2));
E = 2*S*J2*(cos(pi*(h+k))+cos(pi*(h-k)));
F = 1/2*S*JK*xi.*eta;

delta = abs((B-C).*conj(D-F) - conj(B-C).*(D-F));

sum1a = A.^2 - E.^2 + abs(B-C).^2 - abs(D-F).^2;
sum1b = A.^2 - E.^2 + abs(B+C).^2 - abs(D+F).^2;

sum2a = sqrt(4*abs(A.*(B-C)+E.*(D-F)).^2-delta.^2);
sum2b = sqrt(4*abs(A.*(B+C)-E.*(D+F)).^2-delta.^2);

omega = [sqrt(sum1a + sum2a); sqrt(sum1a - sum2a); sqrt(sum1b + sum2b); sqrt(sum1b - sum2b)];

end