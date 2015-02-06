%% Description
% This tutorial reproduces the calculated spin wave spectrum of
% YbLATEX_2PATEXTiLATEX_2PATEXOLATEX_7PATEX with the magnetic Hamiltonian
% proposed in the following paper: <http://journals.aps.org/prx/abstract/10.1103/PhysRevX.1.021002 PRX *1* , 021002 (2011)>. 


%% Create crystal structure
% To create the cubic crystal structure of YbLATEX_2PATEXTiLATEX_2PATEXOLATEX_7PATEX,
% we need to the
% exact lattice parameter is unimportant for the spin wave calculation as
% long as we are using lattice units. The spin of the magnetic atoms are
% automatically created from the ion label that contains the ionic charge
% after the element label. We also define the non-magnetic atoms for
% plotting.

sw_addsym('-z, y+3/4, x+3/4; z+3/4, -y, x+3/4; z+3/4, y+3/4, -x; y+3/4, x+3/4, -z; x+3/4, -z, y+3/4; -z, x+3/4, y+3/4','F d -3 m Z');

ybti = sw;
a = 10.0307;
ybti.genlattice('lat_const',[a a a],'angled',[90 90 90],'sym','F d -3')
ybti.addatom('label','Yb3+','r',[1/2 1/2 1/2])
ybti.addatom('label','Ti4+','r',[0 0 0])
ybti.addatom('label','O2-','r',[0.3318 1/8 1/8])
ybti.addatom('label','O2-','r',[3/8 3/8 3/8])
plot(ybti)

%% Plot cubic environment of Yb

% draw oxygen polyhedra
sw_drawpoly('cAtom',1,'pAtom',3:4,'limits',8);


%% create spin Hamiltonian

% remove non-magnetic atoms, not necessary
ybti.unit_cell.r = ybti.unit_cell.r(:,1);
ybti.unit_cell.S = ybti.unit_cell.S(1);
ybti.unit_cell.label = ybti.unit_cell.label(1);
ybti.unit_cell.color = ybti.unit_cell.color(:,1);

% generate bonds
ybti.gencoupling

% create 3x3 matrices
ybti.addmatrix('label','J1')
ybti.addmatrix('label','g0','value',-0.84*ones(3)+4.32*eye(3));

% assigne J1 to 1st neighbour bonds
ybti.addcoupling('J1',1)
% assigne g0 as g-tensor for every atom
ybti.addg('g0')


% define 2 different magnetic field direction and field strength
n = [1 -1 0];
% 2 different field strength
B1 = 5; % Tesla
B2 = 2; % Tesla

% the anisotropy matrix has the form: [A B B;B A B;B B A]
% 3-fold rotation symmetry along the (111) direction
% its eigenvalues are: g_xy=A-B; g_z = A + 2*B
ybti.matrix.mat(:,:,2) =  -0.84*ones(3)+4.32*eye(3);

% J-values
J1 = -0.09; J2 = -0.22; J3 = -0.29; J4 = 0.01;
% symmetry analysis on the allowed exchange matrix elemnts:
ybti.getmatrix('label','J1');
% assign J1...J4 to the right matrix elements
% note the minus sign in from of J4
ybti.setmatrix('label','J1','pref',[J1 J3 J2 -J4]);

%% calculate spin wave spectrum

% define list of Q-scans
Q = {};
Q{1} = {[-0.5 -0.5 -0.5] [2 2 2]};
Q{2} = {[1 1 -2] [1 1 1.5]};
Q{3} = {[2 2 -2] [2 2 1.5]};
Q{4} = {[-0.5 -0.5 0] [2.5 2.5 0]};
Q{5} = {[0 0 1] [2.3 2.3 1]};

% new ifugre
figure

% loop over list of Q-scans
for ii = 1:10
    
    % select appropriate field direction
    if ii <= 5
        B = B1;
    else
        B = B2;
    end
    % assign magnetic field
    ybti.field(n/norm(n)*B);
    
    if (ii == 1) || (ii==6)
        % create fully polarised magnetic structure along the field direction
        ybti.genmagstr('S',n','mode','helical');
        % find best structure using steepest descendend
        ybti.optmagsteep('nRun',100)
    end
    
    % spin wave spectrum
    ybtiSpec = ybti.spinwave([Q{mod(ii-1,5)+1} {200}],'gtensor',true);
    % neutron scattering cross section
    ybtiSpec = sw_neutron(ybtiSpec);
    % bin the spectrum in energy
    ybtiSpec = sw_egrid(ybtiSpec,'Evect',linspace(0,2,500),'component','Sperp');
    
    % subplot
    subplot(2,5,ii)
    % colorplot with finite energy resolution FWHM 0.09 meV
    sw_plotspec(ybtiSpec,'axLim',[0 0.5],'mode',3,'dE',0.09,'colorbar',false,'legend',false);
    title('')
    % plot dispersion
    %sw_plotspec(ybtiSpec,'axLim',[0 2],'mode',1,'dE',0.09,'colormap',[0 0 0],'colorbar',false);
    %
    caxis([0 60])
    colormap(jet)
end