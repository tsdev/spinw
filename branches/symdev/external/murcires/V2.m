%% V2 parameters for alpha-CaCr2O4 experiment

% small number
delta = 1e-5;
% big number
bignum = 1e5;
% neutron energy transfer in meV
dE = 0;
% fixed neutron wavenumber in Angstrom^-1
kfix = 1.55;
% fixed beam: 1 fix ki, 2 fix kf
fx = 2;
%
switch fx
    case 1
        ki = kfix;
        Ei = sw_converter(ki,'k','meV'); % in meV
        Ef = Ei - dE;
        kf = sw_converter(Ef,'meV','k'); % in Angstrom^-1
    case 2
        kf = kfix;
        Ef = sw_converter(kf,'k','meV'); % in meV
        Ei = Ef + dE;
        ki = sw_converter(Ei,'meV','k'); % in Angstrom^-1
end

% momentum transfer in r.l.u.
qrlu = [1 -1.33 0];
% create crystal lattice
lat = sw;
lat.genlattice('lat_const',[11.057 5.794 5.118],'angled',[90 90 90])
% momentum transfer in Angstrom^-1
qA = qrlu*2*pi*inv(lat.basisvector);
% 2 Theta angle in radian
thsam = acos((ki^2+kf^2-norm(qA)^2)/(2*ki*kf))/2;
% neutron wavelength
lambdai = 2*pi/ki;

% GUIDE
% width in m (from cm)
pars.g = [3 12.5]/100;
% divergence  in std dev min (from min/Angstrom)
% sqrt(12) is the standard deviation of the box function
guidiv = 2*7.05*lambdai*sqrt(12);
% guide monochromator distance in m (from cm)
pars.dgm = 15/100;

% SOLLER COLLIMATORS
% divergence in std. dev. min
pars.sollers = [50 bignum;60 bignum;bignum bignum;bignum bignum];
% use the guide divergence and the soller divergence and take the minimum
% for the first soller
pars.sollers(1,1) = min(pars.sollers(1,1),guidiv);
pars.sollers(1,2) = min(pars.sollers(1,2),guidiv);
% convert divergence to radian (from std. dev. min)
pars.sollers = pars.sollers*pi/180/60;

% SLITS
% position, width
pars.slits = [];

% MONOCHROMATOR
% number of blades
pars.blades = [12 12];
% thickness, width, height in m (from cm)
pars.mon = [0.2 [12.5 12.5]./pars.blades]/100;
% d-spacing monochromator in Angstrom
pars.dmon = 3.355;
% theta monochromator
thmon = asin(pi/pars.dmon/ki);
% focus point 1/r*sin(theta) in m (from 1/m)
pars.foc = [1/delta*sin(thmon) 1/delta*sin(thmon) 1/0.6];
% mosaicity of monochromator crystal in radian FWHM (from min FWHM???)
pars.mos = 30*pi/180/60;
% scattering sense
pars.rots = [1 -1 1];
% incident neutron energy
pars.Ei = Ei;
% monochromator sample distance in m (from cm)
pars.dms = 165/100;

% SAMPLE
% sample shape
pars.sampsh = 'ccc';
% sample dimensions in m (average in the plane) (from cm)
pars.samp = [0.45 0.45 0.5]/100;
% sample shape orientation, xyz principal axes in row vectors
pars.Usamp = eye(3);
% 2 Theta angle in degree
pars.th2 = thsam;
% sample analyser distance in m (from cm)
pars.dsa = 130/100;

% ANALYSER
% number of blades
pars.blades(3:4) = [12 12];
% thickness, width, height in m (from cm)
pars.ana = [0.2 [20 12.5]./pars.blades(3:4)]/100;
% d-spacing analyser in Angstrom
pars.dana = 3.355;
% theta analyser
thana = asin(pi/pars.dana/kf);
% focus point 1/r*sin(theta) in m (from 1/m)
%pars.foc(4:6) = [1/1*sin(thana) 1/1*sin(thana) 1/delta];
% optimised focus
pars.foc(4:6) = [1/0.53442*sin(thana) 1/0.53442*sin(thana) 1/delta];
%pars.foc(4:6) = [1/0.46471*sin(thana) 1/0.46471*sin(thana) 1/delta];
% mosaicity of analyser crystal in radian FWHM (from min FWHM???)
pars.mos(2) = 30*pi/180/60;
% incident neutron energy
pars.Ef = Ef;
% analyser detector distance in m (from cm)
pars.dad = 100/100;

% DETECTOR
% detector shape
pars.detsh = 'bbb';
% detector dimensions in m (from cm)
pars.det = [1 4.9 15]/100;
% detector shape orientation, xyz principal axes in row vectors
pars.Udet = eye(3);

pars = rmfield(pars,'sollers');

pars1 = pars;
pars1.blades = pars1.blades(1:2);
pars1.foc = pars1.foc(1:3);
pars1.mos = pars1.mos(1);

pars2 = pars;
pars2.blades = pars2.blades(3:4);
pars2.foc = pars2.foc(4:6);
pars2.mos = pars2.mos(2);

resMat = fullres(pars1,pars2,'mon','ana',norm(qA),dE,kfix,fx);

%% resolution ellipsoid in the (qx,qy,qz,E) coordinate system
% qx || q, 
% integrate out the energy resolution














