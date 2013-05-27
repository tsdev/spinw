function [Res, errMsg] = sw_test_sw2(tol)

Res = 0;
errMsg = [];

try
    %% LiNiPO4
    % Model taken from T. Jensen et al., Phys. Rev. B. 79, 092413(2009)
    
    
    % Fitted parameters at T = 1.5 K
    
    % Jbc = 1.04;
    % Jb  = 0.670;
    % Jc  = -0.05;
    % Jac = -0.11;
    % Jab = 0.30;
    % Da  = 0.339;
    % Db  = 1.82;
    % Dc  = 0;
    Jbc = 1.036;
    Jb  = 0.6701;
    Jc  = -0.0469;
    Jac = -0.1121;
    Jab = 0.2977;
    Da  = 0.1969;
    Db  = 0.9097;
    Dc  = 0;
    
    Jvect = [Jbc Jab Jb Jc Jac Da Db];
    %% Analytical dispersion from the paper
    
    nQ = 500;
    hFig = figure;
    
    subplot(2,2,1)
    Qa = [0 1 0];
    Qb = [2 1 0];
    omega = Jensen_LiNiPO4({Qa Qb nQ}, Jvect);
    plot(linspace(0,2,nQ),omega(1,:),'b-',linspace(0,2,nQ),abs(omega(2,:)),'r-.');
    axis([0 2 0 8.5]);
    xlabel('H (r.l.u)'); ylabel('\omega (meV)'); title('(H,1,0)');
    
    subplot(2,2,2)
    Qa = [0 0 0];
    Qb = [0 2 0];
    omega = Jensen_LiNiPO4({Qa Qb nQ}, Jvect);
    plot(linspace(0,2,nQ),omega(1,:),'b-',linspace(0,2,nQ),abs(omega(2,:)),'r-.');
    axis([0 2 0 8.5]);
    xlabel('K (r.l.u)'); ylabel('\omega (meV)'); title('(0,K,0)');
    
    subplot(2,2,3)
    Qa = [0 1 0];
    Qb = [0 1 2];
    omega = Jensen_LiNiPO4({Qa Qb nQ}, Jvect);
    plot(linspace(0,2,nQ),omega(1,:),'b-',linspace(0,2,nQ),abs(omega(2,:)),'r-.');
    axis([0 2 0 8.5]);
    xlabel('L (r.l.u)'); ylabel('\omega (meV)'); title('(0,1,L)');
    close(hFig);
    
    %% LiNiPO4 structure in SpinW (only magnetic atoms)
    
    linipo = sw;
    linipo.genlattice('lat_const',[10.02 5.86 4.68],'sym','P n m a')
    
    linipo.addatom('r',[1/4 1/4 0],'S',1,'label','MNi2','color',[0;0;255])
    %Ni2.r = [0.2756; 1/4; 0.9825];
    
    linipo.gencoupling
    
    % Define the interactions
    
    % Isotropic exchanges
    linipo.addmatrix('mat',eye(3)*Jbc,'color',[255 0 0],'label','Jbc')
    linipo.addmatrix('mat',eye(3)*Jb, 'color',[0 255 0],'label','Jb' )
    linipo.addmatrix('mat',eye(3)*Jc, 'color',[0 0 255],'label','Jc' )
    linipo.addmatrix('mat',eye(3)*Jab,'color',[0 125 125],'label','Jab' )
    linipo.addmatrix('mat',eye(3)*Jac,'color',[125 125 0],'label','Jac' )
    
    % anisotropy matrix
    linipo.addmatrix('mat',diag([Da Db Dc]), 'color',[125 0 125],'label','D' )
    
    
    linipo.addcoupling('Jbc',1)
    linipo.addcoupling('Jc' ,2)
    linipo.addcoupling('Jab',[5 6])
    linipo.addcoupling('Jac',[3 4])
    linipo.addcoupling('Jb' ,7)
    
    linipo.addaniso('D')
    
    linipo.genmagstr('mode','direct','S',[0 0 0 0;0 0 0 0;1 -1 -1 1]);
    
    %% plot structure
    
    hFig = plot(linipo,'range',[-0.3 0.8;-0.3 0.8;-0.1 1.1]);
    close(hFig);
    
    %% simulated annealing
    
    %     par_anneal.initT = 100;
    %     par_anneal.endT  = 1e-2;
    %     par_anneal.nMC   = 1000;
    %     par_anneal.cool  = @(T)0.8*T;
    %     par_anneal.nStat = 1;
    %
    %     [linipo, aStat] = sw_anneal(linipo,par_anneal);
    
    
    %% spin wave dispersion
    
    nQ = 400;
    Qa = [0 1 0];
    Qb = [2 1 0];
    
    specLi = linipo.spinwave({Qa Qb nQ});
    specLi = sw_neutron(specLi,'pol',false);
    specLi = sw_conv(specLi,'Evect',linspace(0,8.5,400));
    
    % calculate difference between exact and numerical solutions
    specLiSim = specLi;
    specLiSim.omega = repmat(Jensen_LiNiPO4({Qa Qb nQ},Jvect),4,1);
    specLi    = sw_omegasum(specLi,'zeroint',1e-8);
    specLiSim = sw_omegasum(specLiSim,'zeroint',1e-8);
    omegaDiff = abs(specLi.omega-specLiSim.omega);
    omegaDiff(isnan(omegaDiff)) = 0;
    omegaDiff = max(max(omegaDiff));
    ratioLi1  = omegaDiff/max(max(abs(real(specLi.omega))));
    
    hFig = figure;
    subplot(2,1,1)
    sw_plotspec(specLi,'axLim',1,'mode',4,'aHandle',gca)
    subplot(2,1,2)
    sw_plotspec(specLi,'axLim',15,'mode',2,'aHandle',gca,'imag',false)
    close(hFig)
    
    Qa = [0 0 0];
    Qb = [0 2 0];
    
    specLi = linipo.spinwave({Qa Qb nQ});
    specLi = sw_neutron(specLi,'pol',false);
    specLi = sw_conv(specLi,'Evect',linspace(0,8.5,400));
    
    hFig = figure;
    subplot(2,1,1)
    sw_plotspec(specLi,'axLim',1,'mode',4,'aHandle',gca)
    subplot(2,1,2)
    sw_plotspec(specLi,'axLim',15,'mode',2,'aHandle',gca,'imag',false)
    close(hFig)
    
    Qa = [0 1 0];
    Qb = [0 1 2];
    
    specLi = linipo.spinwave({Qa Qb nQ});
    specLi = sw_neutron(specLi,'pol',false);
    specLi = sw_conv(specLi,'Evect',linspace(0,8.5,400));
    
    % calculate difference between exact and numerical solutions
    specLiSim = specLi;
    specLiSim.omega = repmat(Jensen_LiNiPO4({Qa Qb nQ},Jvect),4,1);
    specLi    = sw_omegasum(specLi,'zeroint',1e-8);
    specLiSim = sw_omegasum(specLiSim,'zeroint',1e-8);
    omegaDiff = abs(specLi.omega-specLiSim.omega);
    omegaDiff(isnan(omegaDiff)) = 0;
    omegaDiff = max(max(omegaDiff));
    ratioLi3  = omegaDiff/max(max(abs(real(specLi.omega))));
    
    hFig = figure;
    subplot(2,1,1)
    sw_plotspec(specLi,'axLim',1,'mode',4,'aHandle',gca)
    subplot(2,1,2)
    sw_plotspec(specLi,'axLim',15,'mode',2,'aHandle',gca,'imag',false)
    close(hFig)
    
catch errMsg
    Res = 1;
    return;
end

try
    if max(abs(ratioLi1),abs(ratioLi3)) > tol
        error('sw_test_sw1:DataError','The calculated spin wave dispersion differs from the exact result!');
    end
catch errMsg
    Res = 2;
    return;
end

end

function omega = Jensen_LiNiPO4(Q, Jinp)
% [omega(1,:)out omega(2,:)out] = Jensen_LiNiPO4(Q, Jinp)
% simulate neutron scattering cross section for CuSb2O6 to compare with sw
%
% Q     Momentum transfer, dimensions are [3 nQ].
% Jinp  Input parameters of the Hamiltonian: [Jyz Jxy Jy Jz Jxz Dx Dy].
%

if iscell(Q)
    Q = sw_qscan(Q);
end

h = Q(1,:);
k = Q(2,:);
l = Q(3,:);

Jyz = Jinp(1);
Jxy = Jinp(2);
Jy  = Jinp(3);
Jz  = Jinp(4);
Jxz = Jinp(5);
Dx  = Jinp(6);
Dy  = Jinp(7);

% spin of the magnetic Ni3+ ions
S  = 1;

A = 4*S*(Jyz + Jxy)-2*S*(Jy*(1-cos(2*pi*k)) + Jz*(1-cos(2*pi*l))+...
    Jxz*(2-cos(pi*(h+l))-cos(pi*(h-l))))+ Dx * S + Dy * S;
B = S*(Dx - Dy);
D = 2*Jyz*S*(cos(pi*(k+l))+cos(pi*(k-l)))+2*Jxy*S*(cos(pi*(h+k))+cos(pi*(h-k)));

omega = [sqrt(A.^2-(B+D).^2); -sqrt(A.^2-(B-D).^2)];

end