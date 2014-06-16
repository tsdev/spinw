function out=Rowlanddata(l1,l2,th)

gamma=2*th;
gammap=sign(th)*180-2*th;
l3=(l1.^2+l2.^2-2*l1.*l2.*cosd(sign(th)*180-gamma)).^0.5;

alpha=asind(l1*sind(gamma)/l3);
beta=2*th-alpha;

ph1=abs(alpha);
ph2=180-abs(beta);

out.alpha=alpha;
out.beta=beta;
out.gamma=gamma;
out.l1=l1;
out.l2=l2;
out.l3=l3;
out.R=l1/sind(abs(alpha))/2;

out.Rm=[l1 0];
out.Rs=[l1+l2*cosd(2*th),l2*sind(2*th)];
out.Rc=[l1/2 l2*sind(2*th)/2+(-l1-l2*cosd(2*th))/2/(-sind(2*th))*cosd(2*th)];
out.Rf=out.Rc+sign(2*th)*out.R*[1 0];
vphm=out.Rm-out.Rc;
vm1=[-vphm(2),vphm(1)]*sign(vphm(1));
out.ph1=atan2(vm1(2),vm1(1))/pi*180;
out.vm1=[cosd(out.ph1),sind(out.ph1)];

out.ph2=out.ph1-2*th;
out.vm2=[cosd(out.ph2),sind(out.ph2)];
out.qph=th+sign(th)*90;
out.Rph=atan2(out.Rc(2)-out.Rm(2),out.Rc(1)-out.Rm(1))/pi*180;

end