function out = monkires(g,dgm,mon,blades,foc,mos,dmon,dms,samp,sampsh,Usamp,Ei,rots,sollers,slits)
% Calculates the resolution matrix for primary instrument with monochromator.
%
% out = MONKIRES(g,dgm,mon,blades,foc,mos,dmon,dms,ord,samp,Ei,rots,sollers, slits)
%
% Input parameters:
%
% Guide:
%
% g         horizontal and vertical dimensions of the quide [y,z] in meter
% dgm       guide end point-monochromator distance in meter
%
% Monochromator:
%
% mon       monochromator crystal paramters: [thickness, width, height] in
%           meter
% blades    number of blades of the monochromator along y and z axis 
%           [N_y N_z], if no focusing then set it to 1.
% foc       monochromator radiuses [f1 f2 r1] in meter:
%               f1  horizontal focal length towards the end point of the guide
%               f2  horizontal focal length towards the sample
%               r1  vertical radius of curvature
%           for flat monochromator, use Inf
% mos       mosaicity of monochromator crystal in radian (FWHM)
% dmon      d-spacing of monochromator in Angstrom
% rots      scattering sense of the monochromator:
%                1  counterclockwise,
%               -1  clockwise
% Ei        incident neutron energy in meV
% dms       monochromator-sample distance in meter
%
% Sample:
%
% sampsh    sample shape along three orthogonal directions, three letters,
%           e.g. 'ccc', for the meaning of each letter see optshape()
%           function help
% samp      sample dimensions along the three orthogonal directions in
%           meter (in most of the cases the principal axis directions of
%           the ellipsoid).
% Usamp     orientation of the sample, in the xyz coordinate system, where
%               x   parallel to the monochromator - sample line
%               y   perpendicular to x in the horizontal plane 
%                   (guide-mono-sample plane)
%               z   vertical
%           matrix with dimensions of [3,3], [v_x; v_y; v_z] three row
%           vectors defining the directions of the three principal axes.
%
% Soller collimator:
%
% sollers   soller collimator parameters [s1 s2; s3 s4] in radian (standard
%           deviation)
%               s1 horizontal before monochromator
%               s2 vertical before monochromator
%               s3 horizontal after monochromator
%               s4 vertical after monochromator
%           Default is Inf(2).
%
% Slits:
%
% slits     slit parameters, matrix with dimensions of [4 N_slit], every
%           column has the following numbers:
%               ind index of position, 1: guide-mono, 2: mono-sample
%               pos position from the previous element in meter
%               w_H horizontal width in meter
%               w_V vertical width in meter
%
% xyz coordinate system for every element is generally the following
% x         parallel to the neutron beam direction
% y         horizontal, perpendicular to x
% z         vertical
%
% Input can be also a single struct type variable with the above fields.
%
% Output parameters:
%
% out is struct type with the following fields:
% pars      matrix with dimensions of [10 13], with the following column 
%           vectors: 
%           [gy,gz,monx1,mony1,monz1,monx2,mony2,monz2,sampx,sampy,sampz,fy,fz]
% sigs      vector of standard deviations in the coordinate system of the 
%           element: [gx gy monx mony monz bladey bladez sampx sampy sampz]
% ds        given widths of the elements with the same order as sigs
% Ri        resolution matrix of the front end in the parameter space,
%           with dimensions of [10 10], if no blade [8 8]
% Ri0       diagonal matrix, equal to diag(sigs)^(-2)
% Rid       other elements of the resolution matrix
% Ci        covariance matrix of parameters, equal to the inverse of Ri
% Cik       covariance matrix in [ki,Ei] space, with dimensions of [4 4]
% Tri       transformation matrix from parameter space to [ki,Ei], with
%           dimensions of [10 4]
% Di        FWHM-s of ki and Ei, effectively the Vanadium width, dimensions
%           of [1 4]
% sizes     real and effective sizes of elements [s_real; s_effective] in
%           meter in the same order as sigs
% dti       delta t, time resolution
% D         soller collimator effect
% S         slit effect
% M         mosaicity effect
%           
%
% The resolution matrix in the parameter space (Rd) has the following
% dimensions ({...} optional):
%
%  1  2   3    4    5  6         7         8       9       10
% gy gz monx mony monz {bladesy} {bladesz} samplex sampley samplez
%
% See also OPTSHAPE, FULLRES, ANAKFRES.
%

if nargin==1
    a=g;
    g=a.g;
    dgm=a.dgm;
    mon=a.mon;
    blades=a.blades;
    foc=a.foc;
    mos=a.mos;
    dmon=a.dmon;
    dms=a.dms;
    samp=a.samp;
    sampsh=a.sampsh;
    if isfield(a,'Usamp')
        Usamp=a.Usamp;
    else
        Usamp=eye(3);
    end;
    Ei=a.Ei;
    rots=a.rots;
    if isfield(a,'sollers')
        sollers=a.sollers;
    else
        sollers=zeros(2);
    end;
    if isfield(a,'slits')
        slits=a.slits;
    else
        slits=0;
    end;
end;

% setting the focuses
if max(blades)==1
    foc=ones(1,3)*10000; % ???
end

if size(foc,2)==1
    foc=ones(1,3)*foc;
end

if size(foc,2)==2
    foc=[foc 1e4];
end

% basic calculations k, lambda mon. th...
ki=NDC('E',Ei,'k');
% neutron wavelength
wli=2*pi/ki;
out.li=dgm+dms;
out.ki=ki;
out.Ei=Ei;
out.vi=NDC('k',ki,'v');
out.ti=out.li/out.vi;

out.sh(1)=optshape(g(1),'b');
out.sh(2)=optshape(g(2),'b');
out.sh(3)=optshape(mon(1),'b');
out.sh(4)=optshape(mon(2),'b');
out.sh(5)=optshape(mon(3),'b');
np=5;
sinmth=wli/2/dmon;
mth=rots(1)*asind(sinmth);

% if no focusing, blades=1
blades=[1 1].*blades;

% checking horizontal and vertical focusing
% kill shows which line is not important if there is no focusing
kill=0;
dang=zeros(1,2);
%dblades=mon(2:3);
if blades(1)<2
    kill=6;
    dang(1)=0;
else
    dang(1)=sinmth*mon(2)/min(foc(1:2));
    np=np+1;
    out.sh(np)=optshape(blades(1),'d');
end;
if blades(2)<2
    kill(sign(kill)+1)=7;
    dang(2)=0;
else
    dang(2)=mon(3)/foc(3);
    np=np+1;
    out.sh(np)=optshape(blades(2),'d');
end;
np=np+1;
out.sh(np)=optshape(samp(1),sampsh(1));
np=np+1;
out.sh(np)=optshape(samp(2),sampsh(2));
np=np+1;
out.sh(np)=optshape(samp(3),sampsh(3));
sigs=zeros(1,np);
ds=zeros(1,np);
for n=1:np
    sigs(n)=out.sh(n).sig;
    ds(n)=out.sh(n).d;
end;

Rowland=Rowlanddata(foc(1),foc(2),mth);
%dblades(1)=sinmth*mon(2)/min(max([sind(Rowland.ph1),sind(Rowland.ph2)]));
%dblades(2)=mon(3);

% if there are  no sollers, then sollers=0
if max(size(sollers))==1
    sollers=ones(2)*100000;
end;

% setting the base vectors
tmp=zeros(np,1);
gy=tmp;
gy(1)=1;
gz=tmp;
gz(2)=1;

%rot=rots(1);

% transformation from the system of mono. cr. to system of ki0
Trmon=[cosd(90-mth) -sind(90-mth) 0;sind(90-mth) cosd(90-mth) 0;0 0 1];
monx1=tmp;
monx1(3:5)=Trmon(:,1);
if min(kill-6)~=0;
    monx1(6)=Rowland.vm1(1)*mon(2);
end;
mony1=tmp;
mony1(3:5)=Trmon(:,2);  %mono y position looking from the direct beam
if min(abs(kill-6))~=0;
    mony1(6)=Rowland.vm1(2)*mon(2);
end;
monz1=tmp;
monz1(3:5)=Trmon(:,3);   % mono z position
if min(abs(kill-7))~=0;
    monz1(7-sum(sign(kill)))=mon(3);
end;

% transformation from the system of mono. cr. to system of ki
Trmon=[cosd(90+mth) -sind(90+mth) 0;sind(90+mth) cosd(90+mth) 0;0 0 1];
monx2=tmp;
monx2(3:5)=Trmon(:,1);
if min(kill-6)~=0;
    monx2(6)=Rowland.vm2(1)*mon(2);
end;
mony2=tmp;
mony2(3:5)=Trmon(:,2); % mono y position in the system of monochromated beam
if min(kill-6)~=0;
    mony2(6)=Rowland.vm2(2)*mon(2);
end;

monz2=monz1;
% orientations due to focusing
fy=tmp;
if min(abs(kill-6))~=0;
    fy(6)=dang(1);
end;
fz=tmp;
if min(abs(kill-7))~=0;
    fz(7-sum(sign(kill)))=dang(2);
end;

sampinds=(8:10)-sum(sign(kill));
sampx=tmp;
sampx(sampinds)=Usamp(:,1); % sample x position
sampy=tmp;
sampy(sampinds)=Usamp(:,2); % sample y position
sampz=tmp;
sampz(sampinds)=Usamp(:,3); % sample z position

% divergences before and after mono
d1y=(-gy+mony1)/dgm;  % direct beam y div.
d1z=(-gz+monz1)/dgm;    % direct beam z div.
d2y=(-mony2+sampy)/dms; % monochr. beam y div.
d2z=(-monz2+sampz)/dms;   % monochr. beam z div.
dli=monx1+monx2+sampx;

out.d=[d1y,d1z,d2y,d2z];
out.dli=dli;

% delta monochromator theta
dmth=(-d1y+d2y)/2; % rots(1) is scattering sense

% hor. angular difference of nominal plane normal and the scattering vector
my=(d1y+d2y)/2-fy; % horizontal deviation from nom. q_mono monfy1/foc(1) is the rotation of q_mono due to focusing
mz=(d1z-d2z)/2/sinmth-fz;
% vert. angular difference of nominal plane normal and the scattering vector

out.dmth=dmth;
out.m=[my,mz];

% delta [ki, E]
dkix=-ki*dmth/tand(mth); % radial change of ki
dkiy=ki*d2y;          % hor. tang. change of ki
dkiz=ki*d2z;          % vert. tang. change of ki
dEi=dkix/ki*Ei*2;       % change of Ei

% setting the diagonal resolutions (due to sizes)
R1=diag(sigs.^-2);

% effect of mosaicity
M=my*my'/mos^2+mz*mz'/mos^2;
out.M=M;
out.m=[my,mz];

% effect of sollers
S=M*0;
D=M*0;
if sollers(1,1)>0
    D=D+d1y*d1y'/sollers(1,1)^2;
else
    sollers(1,1)=1e10;
end;
if sollers(1,2)>0
    D=D+d1z*d1z'/sollers(1,2)^2;
    % else
    %     sollers(1,2)=1e10;
end;
if sollers(2,1)>0
    D=D+d2y*d2y'/sollers(2,1)^2;
    % else
    %     sollers(2,1)=1e10;
end;
if sollers(2,2)>0
    D=D+d2z*d2z'/sollers(2,2)^2;
    % else
    %     sollers(2,2)=1e10;
end;
out.D=D;

if max(size(slits))>1
    tmpp=[gy,gz,mony2,monz2];
    tmpp(:,:,2)=[mony1,monz1,sampy,sampz];
    tmpd=[dgm,dms];
    for ns=1:size(slits,1)
        tmp=slits(ns,1)-1;
        posinds=tmp*2+1;
        if slits(ns,3)>0
            tmppos=tmpp(:,posinds,1)+slits(ns,2)/tmpd(tmp+1)*(tmpp(:,posinds,2)-tmpp(:,posinds,1));
            S=S+tmppos*tmppos'/slits(ns,3)^2/12;
            out.slits(1:np,ns,1)=tmppos;
        end;
        posinds=tmp*2+2;
        if slits(ns,4)>0
            tmppos=tmpp(:,posinds,1)+slits(ns,2)/tmpd(tmp+1)*(tmpp(:,posinds,2)-tmpp(:,posinds,1));
            S=S+tmppos*tmppos'/slits(ns,4)^2/12;
            out.slits(1:np,ns,2)=tmppos;
        end;
        
    end;
end;
out.S=S;
Rd=M+S+D;

out.dti=(out.dli/out.li-dkix/ki)*out.ti;
% resolution matrix in the parameter space
Rp=R1+Rd;

% transformation matrix from the parameter space to [ki,Ei];
Tr=[dkix,dkiy,dkiz,dEi];

% reducing matrices if there is no focusing
out.pars=[gy,gz,monx1,mony1,monz1,monx2,mony2,monz2,sampx,sampy,sampz,fy,fz];

out.sigs=sigs;
out.ds=ds;
out.Ri0=R1;
out.Rid=Rd;
out.Ri=Rp;
out.Ci=inv(Rp); % covariance of parameters
out.Cik=Tr'*out.Ci*Tr; % covariance of [ki,Ei];
out.Tri=Tr;
out.Di=diag(out.Cik*8*log(2))'.^0.5; % FWHM-s of ki and Ei
out.sizes=[diag(R1)'.^-0.5;diag(out.Ci)'.^0.5].*(ones(2,1)*(ds./sigs));

end