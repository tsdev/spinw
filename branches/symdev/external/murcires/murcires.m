function out = murcires(pars1,pars2,prim,sec,q,E,kifix,kifixn)
% Calculates the resolution ellipsoid
%
% out = MURCIRES(pars1,pars2,prim,sec,q,E,kifix,kifixn)
%
% Input parameters:
%
% pars1     parameters of the front end, see the help of the selected front
%           end function
% pars2     parameters of the back end, see the help of the selected back
%           end functon
% prim      type of the font end instrument:
%               'mon'       monochromator, see monkires() function
%               'TOF'       time of flight, see TOFkires() function
% sec       type of the back end instrument:
%               'ana'       analyser, see anakfres() function
%               'CAMEA'     flat cone analyser, see cameakfres() function
% q         momentum transfer in Angstrom^-1
% E         energy transfer in meV
% kifix     value of the fixed ki/kf in Angstrom^-1
% kifixn    selects the fixed neutron momentum:
%               1   fixed ki
%               2   fixed kf
%
% Output parameters:
%
% out struct type contains all output parameter of the front end and back
% end functions and extra fields:
%
% R         resolution matrix in (q,E) space
% Rp        resolution matrix in parameter space
% C         covariance matrix in (q,E) space
% Cp        covariance matrix in parameter space
% DqE       resolution along (q,E) if all other directions are integrated
%           out (effectively Vanadium width)
%
% See also ANAKFRES, MONKIRES.
%

out.in={pars1,pars2,prim,sec,q,E,kifix,kifixn};
tangs=tasqE2angles(q,E,kifix,kifixn,pars2.rots(2));
ki=tangs.ki;
kf=tangs.kf;
Ei=tangs.Ei;
Ef=tangs.Ef;
th2=tangs.th2;
pars1.Ei=tangs.Ei;
TOFi=0;
TOFf=0;
switch prim
    case 'mon'
        pr=monkires(pars1);
%     case 'TOF'
%         pr=TOFkires(pars1);
%         TOFi=1;
end;
pars2.th2=th2;
pars2.Ef=tangs.Ef;
switch sec
    case 'ana'
        ss=anakfres(pars2);
%     case 'CAMEA'
%         ss=cameakfres(pars2);
end;

s1=size(pr.Ri0,1);
s2=size(ss.Rf0,1);

s=s1+s2-3;

tmptr=zeros(s,4);
out.Tri=tmptr;
out.Tri(1:s1,:)=pr.Tri;
if TOFi==1;
    out.Tri(s1-2:s,4)=out.Tri(s1-2:s,4)+ss.dtf/pr.t*2*Ei;
    out.Tri(:,1)=out.Tri(:,4)/Ei*ki/2;
end;
out.Trf=tmptr;
out.Trf(s1-2:s,:)=ss.Trf;
if TOFf==1;
    out.Trf(1:s1,4)=out.Trf(1:s1,4)+pr.dti/ss.t*2*Ef;
    out.Trf(:,1)=out.Trf(:,4)/Ei*ki/2;
end;
% TODO
%rots(2)=pars1.rots(2);
scattangs=thE2angs(th2,E,tangs.Ef);



out.Tr=out.Tri*scattangs.roti-out.Trf*scattangs.rotf;
out.R0=zeros(s);
out.R0(1:s1,1:s1)=pr.Ri0;
out.R0(s1-2:s,s1-2:s)=ss.Rf0;
out.Rd=zeros(s);
out.Rd(1:s1,1:s1)=out.Rd(1:s1,1:s1)+pr.Rid;
out.Rd(s1-2:s,s1-2:s)=out.Rd(s1-2:s,s1-2:s)+ss.Rfd;
out.D=zeros(s);
out.D(1:s1,1:s1)=out.D(1:s1,1:s1)+pr.D;
out.D(s1-2:s,s1-2:s)=out.D(s1-2:s,s1-2:s)+ss.D;
out.S=zeros(s);
out.S(1:s1,1:s1)=out.S(1:s1,1:s1)+pr.S;
out.S(s1-2:s,s1-2:s)=out.S(s1-2:s,s1-2:s)+ss.S;
out.M=zeros(s);
out.M(1:s1,1:s1)=out.M(1:s1,1:s1)+pr.M;
out.M(s1-2:s,s1-2:s)=out.M(s1-2:s,s1-2:s)+ss.M;
nparp=size(pr.pars);
npars=size(ss.pars);
out.pars(1:s1,1:nparp(2))=pr.pars;
out.pars(s1-2:s,nparp(2)+1:nparp(2)+npars(2))=ss.pars;

out.d=zeros(s,8);
out.d(1:s1,1:4)=pr.d;
out.d(s1-2:s,5:8)=ss.d;
out.m=zeros(s,4);
if TOFi==0
    out.m(1:s1,1:2)=pr.m;
end;
if TOFf==0
    out.m(s1-2:s,3:4)=ss.m;
end;
out.Rp=out.R0+out.Rd;
out.Cp=inv(out.Rp);
out.C=out.Tr'*out.Cp*out.Tr;
out.R=inv(out.C);
% sigs=pr.sigs(
out.sizes=[diag(out.R0)'.^-0.5;diag(out.Cp)'.^0.5];
out.DqE=[diag(out.R)'.^-0.5;diag(out.C)'.^0.5]*(8*log(2))^0.5;

tmp=(8*log(2)*[diag(inv(out.Tri(:,1:3)'*out.Cp*out.Tri(:,1:3)))'.^-1;diag(out.Tri(:,1:3)'*out.Cp*out.Tri(:,1:3))']).^0.5;
tmp(:,4)=tmp(:,1)/ki*2*Ei;
out.DkEi=tmp;
tmp=(8*log(2)*[diag(inv(out.Trf(:,1:3)'*out.Cp*out.Trf(:,1:3)))'.^-1;diag(out.Trf(:,1:3)'*out.Cp*out.Trf(:,1:3))']).^0.5;
tmp(:,4)=tmp(:,1)/kf*2*Ef;
out.DkEf=tmp;
out.prim=pr;
out.sec=ss;

end