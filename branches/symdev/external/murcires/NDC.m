function out=NDC(in,dat,outd)
% function NDC(in,dat,outd) changes the data of the neutron 
% in and outd can be 'wl', 'k', 'E' or 'f'
% datas can be: wl (A), k(1/A), v(m/s), E(meV), f(THz)

% constants
Avogadro=6.0221415e23;
kb=1.3806488e-23;
Tstandard=298.15;
Vmolstandard=1*kb*Tstandard*Avogadro/1e5; % [m^3]

h=6.626068e-34;
hv=h/2/pi;
e=1.60217646e-19;
mn=1.674927351e-27;
clight=299792458;
EmeV2J=e/1000;
EJ2meV=1000/e;
EJ2T=1/kb;
ET2J=kb;
EmeV2T=EmeV2J*EJ2T;
ET2meV=ET2J*EJ2meV;
dGr002=6.7079/2;
dGr100=2.4612;
EmeV2om=EmeV2J/hv*1e-12;
Eom2meV=1/EmeV2om;
EmeV2f=EmeV2om/2/pi;
Ef2meV=1/EmeV2f;
EmeVm12ns=1/EmeV2f/1000;
Ensm12meV=1/EmeVm12ns;
EmeV2cmm1=EmeV2f/clight/100*1e12;
Ecmm12meV=1/EmeV2cmm1;

pkbar2GPa=0.1;
pGPa2kbar=10;

sigmaHe=5330*1e-28;
RTmeV=25.8520;
RTwl=1.7789;
SigmaHe1bar=sigmaHe*Avogadro/Vmolstandard; % [1/m] for 25meV neutrons

dGr002=3.3540;
% calculations

    switch in
        case 'wl'
            k=dat.^-1*2*pi;
        case 'k'
            k=dat;
        case 'v'
            k=dat*mn/hv/1e10;
        case 'E'
            k=(dat*EmeV2J*2*mn).^0.5/hv/1e10;
        case 'f'
            k=(dat*Ef2meV*EmeV2J*2*mn).^0.5/hv/1e10;
        case 'PGth'
            k=2*pi./(2*dGr002*sind(dat));
    end;
    
    switch outd
        case 'wl'
            out=k.^-1*2*pi;
        case 'k'
            out=k;
        case 'v'
            out=hv*k/mn*1e10;
        case 'E'
            out=k.^2*hv^2/2/mn*EJ2meV*1e20;
        case 'f'
            out=k.^2*hv^2/2/mn*EJ2meV*EmeV2f*1e20;
        case 'PGth'
            out=asind(2*pi./k/2/dGr002);
    end;
            
end
