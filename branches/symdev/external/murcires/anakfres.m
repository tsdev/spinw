function out=anakfres(samp,sampsh,dsa,ana,blades,foc,mos,dana,dad,det,detsh,Usamp,Udet,Ef,th2,rots,sollers,slits)
% Calculates the resolution matrix for primary instrument with monochromator.
%
% out = ANAKFRES(samp,sampsh,dsa,ana,blades,foc,mos,dana,dad,det,detsh,Usamp,Udet,Ef,th2,rots,sollers,slits)
%
% Input parameters:
%
% Sample:
%
% sampsh    sample shape along three orthogonal directions, three letters,
%           e.g. 'ccc', for the meaning of each letter see optshape()
%           function help
% Usamp     orientation of the sample, in the xyz coordinate system, where
%               x   parallel to the analyser - sample line
%               y   perpendicular to x in the horizontal plane 
%                   (guide-mono-sample plane)
%               z   vertical
%           matrix with dimensions of [3,3], [v_x; v_y; v_z] three row
%           vectors defining the directions of the three principal axes.
% samp      sample dimensions along the three orthogonal directions in
%           meter (in most of the cases the principal axis directions of
%           the ellipsoid).
% dsa       sample - analyser distance in meter
%
% Analyser:
%
% ana       analyser crystal paramters: [thickness, width, height] in
%           meter
% blades    number of blades of the analyser along y and z axis [N_y N_z],
%           if no focusing then set it to 1.
% foc       analyser radiuses [f1 f2 r1] in meter:
%               f1  horizontal focal length towards the sample
%               f2  horizontal focal length towards the detector
%               r1  vertical radius of curvature
%           for flat analyser, use Inf
% mos       mosaicity of analyser crystal in radian (FWHM)
% dana      d-spacing of analyser in Angstrom
% dad       analyser-detector distance in meter
%
% Detector:
%
% detsh    detector shape along three orthogonal directions, three letters,
%           e.g. 'ccc', for the meaning of each letter see optshape()
%           function help
% det       detector dimensions along the three orthogonal directions in
%           meter (in most of the cases the principal axis directions of
%           the ellipsoid).
% Udet     orientation of the detector, in the xyz coordinate system, where
%               x   parallel to the analyser - detector line
%               y   perpendicular to x in the horizontal plane 
%                   (sample-analyser-detector plane)
%               z   vertical
%           matrix with dimensions of [3,3], [v_x; v_y; v_z] three row
%           vectors defining the directions of the three principal axes.
%
% Ef        final neutron energy in meV
% th2       two theta angle in degree (A4 angle)
% rots      scattering sense at the monochromator, sample and analyser 
%           [rm rs ra] with the following values:
%                1  counterclockwise,
%               -1  clockwise
% sollers   soller collimator parameters [s1 s2; s3 s4] in radian (standard
%           deviation)
%               s1 horizontal before analyser
%               s2 vertical before analyser
%               s3 horizontal after analyser
%               s4 vertical after analyser
%           Default is Inf(2).
% slits     slit parameters, matrix with dimensions of [4 N_slit], every
%           column has the following numbers:
%               ind index of position, 1: guide-mono, 2: mono-sample, 3:
%                   sample-ana, 4: ana-detector
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
% Rf        resolution matrix of the back end in the parameter space,
%           with dimensions of [10 10], if no blade [8 8]
% Rf0       diagonal matrix, equal to diag(sigs)^(-2)
% Rfd       other elements of the resolution matrix
% Cf        covariance matrix of parameters, equal to the inverse of Rf
% Cfk       covariance matrix in [kf,Ef] space, with dimensions of [4 4]
% Trf       transformation matrix from parameter space to [kf,Ef], with
%           dimensions of [10 4]
% Df        FWHM-s of kf and Ef, effectively the Vanadium width, dimensions
%           of [1 4]
% sizes     real and effective sizes of elements [s_real; s_effective] in
%           meter in the same order as sigs
% dtf       delta t, time resolution
% D         soller collimator effect
% S         slits effect
% M         mosaicity effect
%           
% The resolution matrix in the parameter space (Rd) has the following
% dimensions ({...} optional):
%
%  1       2       3       4    5    6    7    8    9    10        11
%  samplex sampley samplez monx mony monz detx dety detz {bladesy} {bladesz}
%
% See also OPTSHAPE, FULLRES, MONKIRES.
%


if nargin==1
    a=samp;
    samp=a.samp;
    sampsh=a.sampsh;
    dsa=a.dsa;
    ana=a.ana;
    blades=a.blades;
    foc=a.foc;
    mos=a.mos;
    dana=a.dana;
    dad=a.dad;
    det=a.det;
    detsh=a.detsh;
    if isfield(a,'Usamp')
        Usamp=a.Usamp;
    else
        Usamp=eye(3);
    end;
    if isfield(a,'Udet')
        Udet=a.Udet;
    else
        Udet=eye(3);
    end;
    Ef=a.Ef;
    rots=a.rots;
    th2=abs(a.th2)*rots(2);
    
    if isfield(a,'sollers')
        sollers=a.sollers;
    else
        sollers=zeros(4);
    end;
    if isfield(a,'slits')
        slits=a.slits;
    else
        slits=0;
    end;
end;


out.shapes(1)=optshape(samp(1),sampsh(1));
out.shapes(2)=optshape(samp(2),sampsh(2));
out.shapes(3)=optshape(samp(3),sampsh(3));
out.shapes(4)=optshape(ana(1),'b');
out.shapes(5)=optshape(ana(2),'b');
out.shapes(6)=optshape(ana(3),'b');
out.shapes(7)=optshape(det(1),detsh(1));
out.shapes(8)=optshape(det(2),detsh(2));
out.shapes(9)=optshape(det(3),detsh(3));
np=9;
kill=0;
if blades(1)>1
    np=np+1;
    out.shapes(np)=optshape(blades(1),'d');
else
    kill=np+1;
end;
if blades(2)>1
    np=np+1;
    out.shapes(np)=optshape(blades(2),'d');
else
    kill(sign(kill)+1)=11;
end;

sigs=zeros(1,np);
ds=zeros(1,np);
for nnp=1:np
    sigs(nnp)=out.shapes(nnp).sig;
    ds(nnp)=out.shapes(nnp).d;
end;
out.sigs=sigs;
out.ds=ds;

if max(blades)==1
    foc=ones(1,3)*10000;
end;
if max(size(foc))==1
    foc=foc*ones(1,3);
end;

kf=NDC('E',Ef,'k');
wlf=2*pi/kf;
out.Ef=Ef;
out.kf=kf;
out.lf=dsa+dad;
out.vf=NDC('k',kf,'v');
out.tf=out.lf/out.vf;
sinath=wlf/2/dana;
ath=rots(3)*asind(sinath);

dang(1)=sinath*ana(2)/min(foc(1:2));
dang(2)=ana(3)/foc(3);
Rowland=Rowlanddata(foc(1),foc(2),ath);
%dblades(1)=sinath*ana(2)/min(abs([sind(Rowland.ph1),sind(Rowland.ph2)]));
%dblades(2)=ana(3);

tmp=zeros(np,1);
Rs=[cosd(th2) -sind(th2) 0;sind(th2) cosd(th2) 0; 0 0 1];
Rsamp=Usamp*Rs;
sampx=tmp;
sampx(1:3)=Rsamp(:,1);
sampy=tmp;
sampy(1:3)=Rsamp(:,2);
sampz=tmp;
sampz(1:3)=Rsamp(:,3);


tmpth=(90-ath);
Rana=[cosd(tmpth) -sind(tmpth) 0; sind(tmpth) cosd(tmpth) 0; 0 0 1];
anax1=tmp;
anax1(4:6)=Rana(:,1);
if min(abs(kill-10))~=0;
    anax1(10)=Rowland.vm1(1)*ana(2);
end;
anay1=tmp;
anay1(4:6)=Rana(:,2);
if min(abs(kill-10))~=0;
    anay1(10)=Rowland.vm1(2)*ana(2);
end;
anaz1=tmp;
anaz1(4:6)=Rana(:,3);
if min(abs(kill-11))~=0;
    anaz1(np)=ana(3);
end;


tmpth=(90+ath);
Rana=[cosd(tmpth) -sind(tmpth) 0; sind(tmpth) cosd(tmpth) 0; 0 0 1];
anax2=tmp;
anax2(4:6)=Rana(:,1);
if min(abs(kill-10))~=0;
    anax2(10)=Rowland.vm2(1)*ana(2);
end;
anay2=tmp;
anay2(4:6)=Rana(:,2);
if min(abs(kill-10))~=0;
    anay2(10)=Rowland.vm2(2)*ana(2);
end;
anaz2=tmp;
anaz2(4:6)=Rana(:,3);
if min(abs(kill-11))~=0;
    anaz2(np)=ana(3);
end;

detx=tmp;
detx(7:9)=Udet(:,1);
dety=tmp;
dety(7:9)=Udet(:,2);
detz=tmp;
detz(7:9)=Udet(:,3);

fy=tmp;
if min(abs(kill-10))~=0
    fy(10)=dang(1);
end;
fz=tmp;
if min(abs(kill-11))~=0
    fz(np)=dang(2);
end;

out.pars=[sampx,sampy,sampz,anax1,anay1,anaz1,anax2,anay2,anaz2,detx,dety,detz,fy,fz];
d3y=(-sampy+anay1)/dsa;
d3z=(-sampz+anaz1)/dsa;
d4y=(-anay2+dety)/dad;
d4z=(-anaz2+detz)/dad;

out.d=[d3y,d3z,d4y,d4z];
dlf=sampx+anax1+anax2+detx;
out.dlf=dlf;


dath=(-d3y+d4y)/2;
my=(d3y+d4y)/2-fy;
mz=(d3z-d4z)/2/sinath-fz;
out.m=[my,mz];

dkfx=-kf*dath/tand(ath);
dkfy=kf*d3y;
dkfz=kf*d3z;
dEf=dkfx/kf*Ef*2;

R1=diag(sigs.^-2);
M=my*my'/mos^2+mz*mz'/mos^2;
S=M*0;
D=M*0;
if max(size(sollers))>1
    if sollers(1,1)>0
        D=D+d3y*d3y'/sollers(1,1)^2;
    end;
    if sollers(1,2)>0
        D=D+d3z*d3z'/sollers(1,2)^2;
    end;
    if sollers(2,1)>0
        D=D+d4y*d4y'/sollers(2,1)^2;
    end;
    if sollers(1,2)>0
        D=D+d4z*d4z'/sollers(2,2)^2;
    end;
end;
if max(size(slits))>1
    tmpp=[sampy,sampz,anay,anaz,dety,detz];
    tmpd=[dsa,dad];
    for ns=1:size(slits,1)
        tmp=slits(ns,1)-3;
        posinds=[tmp*2+1,tmp*2+3];
        if slits(ns,3)>0
            tmppos=tmpp(:,posinds(1))+slits(ns,2)/tmpd(tmp+1)*(tmpp(:,posinds(2))-tmpp(:,posinds(1)));
            S=S+tmppos*tmppos'/slits(ns,3)^2;
            out.slits(1:np,ns,1)=tmppos;
        end;
        posinds=[tmp*2+2,tmp*2+4];
        if slits(ns,4)>0
            tmppos=tmpp(:,posinds(1))+slits(ns,2)/tmpd(tmp+1)*(tmpp(:,posinds(2))-tmpp(:,posinds(1)));
            S=S+tmppos*tmppos'/slits(ns,4)^2;
            out.slits(1:np,ns,2)=tmppos;
        end;
        
    end;
end;
Rd=M+S+D;

out.dtf=(dlf/out.lf-dkfx/out.kf)*out.tf;
R=R1+Rd;

Tr=[dkfx,dkfy,dkfz,dEf];
out.D=D;
out.S=S;
out.M=M;
out.Rf0=R1;
out.Rfd=Rd;
out.Rf=R;
out.Cf=inv(R);
out.Cfk=Tr'*out.Cf*Tr;
out.Trf=Tr;
out.Df=diag(out.Cfk*8*log(2)).^0.5;
out.sizes=[diag(inv(R1))';diag(inv(R))'].^0.5.*(ones(2,1)*(ds./sigs));
% out.dl=dls;

end