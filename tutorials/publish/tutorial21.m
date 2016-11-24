%% Description
% This script calculates the low temperature spin dynamics of yttrium iron
% garnet (YIG) to compare with the recently published calculation of J.
% Baker and G. E. W. Bauer,
% <http://journals.aps.org/prl/abstract/10.1103/PhysRevLett.117.217201 PRL *117*, 217201 (2016)>.
% The published calculation is based on a classical spin dynamics described
% by the Landau-Lifshitz-Gilbert equation and a numerical integration in
% time. This method can be used to calculate spin dynamics at finite
% temperature. On the other hand SpinW can calculate the spin wave spectrum
% by diagonalizing the spin Hamiltonian that should be identical to the
% classical solution up to a finite zero point energy shift at low
% temperatures. Here we reproduce the low temperature spectrum to compare
% with the classical results of the paper. The below code will only work
% with the unreleased version of SpinW, that can be found on the gitHub
% repository. 

%% Importing and plotting the YIG crystal
% We get the crystal structure of YIG from a .cif file stored online. SpinW
% is able to download the file from a given link and create the crystal
% structure.

yig = spinw('https://drive.google.com/uc?export=download&id=0BzFs7CQXhehScVFfbkhrZHo1Z1k');

% The imported .cif file contains all symmetry operators of its space
% group. SpinW will determine the generators of the space group and stores
% them in the yig.lattice.sym matrix that has dimensions of 3x4xnOp. Where
% the yig.lattice.sym(:,1:3,i) is the rotation operator of the ith symmetry
% operator, while the yig.lattice.sym(:,4,i) is the corresponding
% translation vector. The imported structure assigns automatically
% different colors for different atoms, thus the same color for the two Fe
% sublattice. Here we color differently the two sublattice:

yig.unit_cell.color(:,3) = sw_colorname('red');
yig.unit_cell.color(:,2) = sw_colorname('blue');

% problem with the generators, the calculated generators from the .cif file
% give only 64 positions, while there should be 96 positions:
% http://it.iucr.org/Ab/ch7o1v0001/sgtable7o1o230/
yig.lattice.sym = sw_gensym('I a -3 d');

% spin quantum number of Fe3+ ions, determined automatically by SpinW
S0 = max(yig.unit_cell.S);

% normalize spins to S=1 as it is in the paper
yig.unit_cell.S = yig.unit_cell.S/S0;

hFig = plot(yig,'pNonMagAtom',0,'labelAtom',0);

% new basis vectors in rows
pBV = [1/2 1/2 -1/2;-1/2 1/2 1/2;1/2 -1/2 1/2];
% lattice constant of YIG
lat = yig.abc(1);

% show the primitive cell basis vectors, plotting unit is Angstrom
hArrow       = sw_arrow([0 0 0],pBV(1,:)*lat,0.06,15,0.5,50);
hArrow(5:8)  = sw_arrow([0 0 0],pBV(2,:)*lat,0.06,15,0.5,50);
hArrow(9:12) = sw_arrow([0 0 0],pBV(3,:)*lat,0.06,15,0.5,50);
set(hArrow,'FaceColor','r')

% add the arrow to the crystal model
sw_addobject(hFig,hArrow);

%% create spin Hamiltonian
% change from BCC to primitive cubic cell
T = yig.newcell({pBV(1,:) pBV(2,:) pBV(3,:)});

% exchange values from the paper
Jad = sw_converter(9.60e-21,'J','THz','photon');
Jdd = sw_converter(3.24e-21,'J','THz','photon');
Jaa = sw_converter(0.92e-21,'J','THz','photon');

% scale the interactions from classical moment size to quantum model
Scl = sqrt(S0*(S0+1));
yig.quickham([Jad Jdd Jaa]/Scl)

% add external field and convert from the standard SpinW unit (meV) to THz
yig.field([0 0 0.01]*sw_converter(1,'meV','THz','photon'))
yig.optmagsteep
yig.genmagstr('mode','rotate','n',[0 0 1])
yig.optmagsteep

if sum(yig.mag_str.F(3,:),2)<0
    yig.mag_str.F = -yig.mag_str.F;
end


%plot(yig,'pnonmagatom',0)


%% Spin wave dispersion to compare with the paper
c
Q0  = T*[1 2 3]';
Q_N = T*[ 1/2  1/2    0]'+Q0;
Q_G = T*[   0    0    0]'+Q0;
Q_H = T*[   0    0    1]'+Q0;

spec = yig.spinwave({Q_N Q_G Q_H 501});

%%

spec = sw_egrid(spec,'component','Sxy-Syx','Evect',linspace(0,28,501));
spec = sw_instrument(spec,'dE',0.75);

figure
sw_plotspec(spec,'mode','disp','colormap',[0 0 0])
hold on
sw_plotspec(spec,'mode','color','imag',true)
colormap(sw_cbrewer('RdBu'))
title('YIG low temperature spin wave spectrum')
ylabel('Energy (THz)')

legend off
colorbar off
set(gca,'XTickLabel',{'N' '\Gamma' 'H'})
caxis([-0.05 0.05])



