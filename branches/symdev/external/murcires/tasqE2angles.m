function out=tasqE2angles(q,E,kfix,fixn,rot)

if max(size(rot))==3
    rot=rot(2);
end;
switch fixn
    case 1
        ki=kfix;
        Ei=NDC('k',ki,'E');
        Ef=Ei-E;
        kf=NDC('E',Ef,'k');
    case 2
        kf=kfix;
        Ef=NDC('k',kf,'E');
        Ei=E+Ef;
        ki=NDC('E',Ei,'k');
    case -1
        Ei=kfix;
        ki=NDC('E',Ei,'k');
        Ef=Ei-E;
        kf=NDC('E',Ef,'k');
    case -2
        Ef=kfix;
        kf=NDC('E',Ef,'k');
        Ei=Ef+E;
        ki=NDC('E',Ei,'k');
end;
if min([Ei,Ef,ki+kf-q])<0
    error('triangle cannot be closed');
end;

% ki=NDC('E',Ei,'k');
% if ki+kf<q
%     error('Cannot close the scattering triangle');
%     return;
% end;
th2=rot*acosd((ki.^2+kf.^2-q.^2)./ki./kf/2);
if abs(real(th2))<1e-5
    if abs(imag(th2))<1e-5
        th2=0;
    else
        error('cannot close the triangle');
    end;
end;


phi=rad2deg(atan2(-kf*sind(th2),ki-kf*cosd(th2)));
psi=phi-th2;
out.ki=ki;
out.Ei=Ei;
out.kf=kf;
out.Ef=Ef;
out.E=E;
out.q=q;
out.th2=th2;
out.phi=phi;
out.psi=psi;
out.roti=[cosd(phi) -sind(phi) 0 0;sind(phi) cosd(phi) 0 0; 0 0 1 0;0 0 0 1];
out.rotf=[cosd(psi) -sind(psi) 0 0;sind(psi) cosd(psi) 0 0; 0 0 1 0;0 0 0 1];

end

function degout = rad2deg(radin)

degout = radin * 180/pi;

end