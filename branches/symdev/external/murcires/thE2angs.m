function out=thE2angs(th2,E,Ef)

kf=NDC('E',Ef,'k');
Ei=Ef+E;
if Ei<0
    error('Ei<0');
    return;
end;

ki=NDC('E',Ei,'k');
q2=ki^2+kf^2-2*ki*kf*cosd(th2);
q=q2^0.5;

out.ki=ki;
out.Ei=Ei;
out.kf=kf;
out.Ef=Ef;
out.E=E;
out.q=q;

s2theta=th2;
out.th2=th2;
phi=atan2((-kf*sind(s2theta)),(ki-kf*cosd(s2theta)));
ph=rad2deg(phi);
out.phi=ph;
roti=[cosd(ph) -sind(ph) 0 0;sind(ph) cosd(ph) 0 0;0 0 1 0;0 0 0 1];
ph=ph-th2;
out.psi=ph;
out.roti=roti;
out.rotf=[cosd(ph) -sind(ph) 0 0;sind(ph) cosd(ph) 0 0;0 0 1 0;0 0 0 1];
out.rotq2k=out.roti';

end

function radout = rad2deg(degin)

radout = degin * pi/180;

end
