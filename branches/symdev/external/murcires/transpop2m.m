function [outm]=transpop2m(f,q,p,mon_flag,nblades,slits)

% conversion from minutes of arc to radians.
pit = pi/180/60; 
% conversion from std. dev. to FWHM of Gaussian
tr=(8*log(2))^0.5;

prim.dmon = p(1);            % monochromator d-spacing in Angs.
sec.dana  = p(2);            % analyser d-spacing in Angs.
prim.mos  = p(3)*pit/tr;     % monochromator mosaic (converted from min FWHM -> radian std. dev.)
% ??? vertical mosaic of the monochromator.
sec.mos   = p(4)*pit/tr;      % analyser mosaic.
% mosana(2)=mosana;     % vertical mosaic spread of the analyser.
% etas=p(5)*pit;        % sample mosaic.
% etasp=etas;	        % vertical mosaic spread of the sample.
rots(1) = p(6);            % scattering sense of monochromator (left=+1,right=-1)
rots(2) = p(7);            % scattering sense of sample (left=+1,right=-1)
rots(3) = p(8);            % scattering sense of analyser (left=+1,right=-1)
prim.slits = slits;
prim.rots  = rots;
sec.rots   = rots;
kfix  = p(9);            % fixed momentum component in ang-1.
kfixn = p(10);           % fx=1 for fixed incident and 2 for scattered wavevector.

soll(1,1) = p(11)*pit/tr;     % ??? horizontal pre-monochromator collimation.
soll(2,1) = p(12)*pit/tr;     % horizontal pre-sample collimation.
soll(3,1) = p(13)*pit/tr;     % horizontal post-sample collimation.
soll(4,1) = p(14)*pit/tr;     % horizontal post-analyser collimation.
soll(1,2) = p(15)*pit/tr;     % vertical pre-monochromator collimation.
soll(2,2) = p(16)*pit/tr;     % vertical pre-sample collimation.
soll(3,2) = p(17)*pit/tr;     % vertical post-sample collimation.
soll(4,2) = p(18)*pit/tr;     % vertical post-analyser collimation.
prim.sollers = soll(1:2,:);
sec.sollers  = soll(3:4,:);
E  = p(34);            % energy transfer.
Es = zeros(1,2);
Es(kfixn) = NDC('k',kfix,'E');
Es(mod(kfixn,2)+1) = sum(Es)+(-1)^kfixn*E;
Ei = Es(1);
Ef = Es(2);
ks = NDC('E',Es,'k');
ki = ks(1);
kf = ks(2);
tmpangs = tasqE2angles(q,E,kfix,kfixn,rots);
sec.th2 = tmpangs.th2;

mth = asind((2*pi/ks(1))/(2*prim.dmon));
ath = asind((2*pi/ks(2))/(2*sec.dana));

% _____________________Extra Parameters________________________________________
offset=42;

sh=['c','b'];
shapesourc=sh(p(1+offset)+1)*[1 1];       % =0 for circular source, =1 for rectangular source.
prim.g(1)=p(2+offset)/100;                % width/diameter of the source (cm).
prim.g(2)=p(3+offset)/100;                % height/diameter of the source (cm).

flag_guide=p(4+offset);          % =0 for no guide, =1 for guide.
guide_h=p(5+offset);             % horizontal guide divergence (mins/Angs)
guide_v=p(6+offset);             % vertical guide divergence (mins/Angs)
if flag_guide
    prim.sollers(1,:)=deg2rad([guide_h,guide_v]/60*2*pi/ks(1))/tr;
end;

tmpsh=sh(p(7+offset)+1);                % =0 for cylindrical sample, =1 for cuboid sample.
prim.sampsh=[tmpsh tmpsh 'b'];
prim.samp(1)=p(8+offset)/100;                % sample width/diameter perp. to Q (cm).
prim.samp(2)=p(9+offset)/100;                % sample width/diameter along Q (cm). 
prim.samp(3)=p(10+offset)/100;               % sample height (cm).

sec.sampsh=prim.sampsh;
sec.samp=prim.samp;


detsh=sh(p(11+offset)+1);               % =0 for circular detector, =1 for rectangular detector.
sec.detsh=[detsh,'b','b'];
detsiz(1)=0.01;
detsiz(2)=p(12+offset)/100;               % width/diameter of the detector (cm).
detsiz(3)=p(13+offset)/100;               % height/diameter of the detector (cm).
sec.det=detsiz;

mon(1)=p(14+offset)/100;               % thickness of monochromator (cm).
% mon(2)=bladesizes(1,1);
% mon(3)=bladesizes(1,2);
% blades(1,1)=p(15+offset)/100/mon(2);
% blades(1,2)=p(16+offset)/100/mon(3);

 mon(2)=p(15+offset)/100/nblades(1,1);               % width of monochromator (cm).
 mon(3)=p(16+offset)/100/nblades(1,2);               % height of monochromator (cm).
prim.mon=mon;
prim.blades=nblades(1,:);

ana(1)=p(17+offset)/100;               % thickness of analyser (cm).


% ana(2)=bladesizes(2,1);
% ana(3)=bladesizes(2,2);
% blades(2,1)=p(18+offset)/100/ana(2);
% blades(2,2)=p(19+offset)/100/ana(3);

ana(2)=p(18+offset)/100/nblades(2,1);               % width of analyser (cm).
ana(3)=p(19+offset)/100/nblades(2,2);               % height of analyser (cm).

sec.ana=ana;
sec.blades=nblades(2,:);




prim.dgm=p(20+offset)/100;                 % distance between source and monochromator (cm).
prim.dms=p(21+offset)/100;                 % distance between monochromator and sample (cm).
sec.dsa=p(22+offset)/100;                 % distance between sample and analyser (cm).
sec.dad=p(23+offset)/100;                 % distance between analyser and detector (cm).

tmpmonfoc=p(24+offset)^-1; % horizontal curvature of monochromator 1/radius (cm-1).
monfoc(1:2)=tmpmonfoc*sind(mth)*[1 1];
monfoc(3)=p(25+offset)^-1;    % vertical curvature of monochromator (cm-1).
prim.foc=monfoc;

tmpanafoc=p(26+offset)^-1; % horizontal curvature of analyser (cm-1).
anafoc(1:2)=tmpanafoc*sind(ath)*[1 1];
anafoc(3)=p(27+offset)^-1;    % vertical curvature of analyser (cm-1).
sec.foc=anafoc;

prim.Ei=Ei;
prim.ki=ki;
sec.Ef=Ef;
sec.kf=kf;

[R0,NP,vi,vf,Error]=rc_popma(f,q,p,mon_flag);
out1=fullres(prim,sec,'mon','ana',q,E,kfix,kfixn);
outm.p=p;
outm.primpar=prim;
outm.secpar=sec;
outm.Rpop=NP;
outm.Rm=out1.R;
outm.full=out1;
outm.pop=R0;

end

function radout = deg2rad(degin)
radout = degin*pi/180;
end