function data=addlinemarker(data,loc_scene,Obj)
a=data.tags.maxsize*Obj.MarkerSize*0.001;
V=Obj.V;
x=V(:,1);
y=V(:,2);
z=V(:,3);

x1=x-a;
x2=x+a;
x3=x-a*0.5;
x4=x+a*0.5;

y1=y-a;
y2=y+a;
y3=y-a*0.5;
y4=y+a*0.5;

z1=z-a;
z2=z+a;


switch(Obj.Marker(1))
    case '*'
        [E,V]=addlinepiece(x,x,y1,y2,z,z);
        [E,V]=addlinepiece(x1,x2,y,y,z,z,E,V);
        [E,V]=addlinepiece(x3,x4,y3,y4,z,z,E,V);
        [E,V]=addlinepiece(x4,x3,y3,y4,z,z,E,V);
    case 'x'
        [E,V]=addlinepiece(x3,x4,y3,y4,z,z);
        [E,V]=addlinepiece(x4,x3,y3,y4,z,z,E,V);
    case '+'
        [E,V]=addlinepiece(x,x,y1,y2,z,z);
        [E,V]=addlinepiece(x1,x2,y,y,z,z,E,V);
        [E,V]=addlinepiece(x,x,y,y,z1,z2,E,V);
    case 's'
        [E,V]=addlinepiece(x1,x1,y1,y2,z,z);
        [E,V]=addlinepiece(x2,x2,y1,y2,z,z,E,V);
        [E,V]=addlinepiece(x1,x2,y1,y1,z,z,E,V);
        [E,V]=addlinepiece(x1,x2,y2,y2,z,z,E,V);
    case 'd'
        [E,V]=addlinepiece(x1,x,y,y1,z,z);
        [E,V]=addlinepiece(x,x2,y1,y,z,z,E,V);
        [E,V]=addlinepiece(x1,x,y,y2,z,z,E,V);
        [E,V]=addlinepiece(x,x2,y2,y,z,z,E,V);
    case 'o',
        [E,V]=addlinepiece(x1,x1,y3,y4,z,z);
        [E,V]=addlinepiece(x2,x2,y3,y4,z,z,E,V);
        [E,V]=addlinepiece(x3,x4,y1,y1,z,z,E,V);
        [E,V]=addlinepiece(x3,x4,y2,y2,z,z,E,V);
        
        [E,V]=addlinepiece(x1,x3,y4,y2,z,z,E,V);
        [E,V]=addlinepiece(x1,x3,y3,y1,z,z,E,V);
        [E,V]=addlinepiece(x2,x4,y3,y1,z,z,E,V);
        [E,V]=addlinepiece(x2,x4,y4,y2,z,z,E,V);
    case '.'
        V= [0.8507 0.5257 0; -0.8507 0.5257 0; -0.8507 -0.5257 0; 0.8507 -0.5257 0;
            0.5257 0 0.8507; 0.5257 0 -0.8507;  -0.5257 0 -0.8507; -0.5257 0 0.8507;
            0 0.8507 0.5257; 0 -0.8507 0.5257;  0 -0.8507 -0.5257; 0 0.8507 -0.5257]*a;
        F=[ 5 8 9; 5 10 8; 6 12 7; 6 7 11; 1 4 5; 1 6 4; 3 2 8; 3 7 2;
            9 12 1; 9 2 12; 10 4 11; 10 11 3; 9 1 5; 12 6 1; 5 4 10; 6 11 4;
            8 2 9; 7 12 2; 8 10 3; 7 3 11];
        F=permute(reshape(bsxfun(@plus,F(:),((0:size(x,1)-1)*size(V,1))),size(F,1),3,[]),[1 3 2]);
        F=reshape(F,size(F,1)*size(F,2),3);
        
        x=bsxfun(@plus,V(:,1),x');
        y=bsxfun(@plus,V(:,2),y');
        z=bsxfun(@plus,V(:,3),z');
        V=[x(:) y(:) z(:)];
    otherwise
        return
end

if(~ischar(Obj.MarkerEdgeColor)), Color=Obj.MarkerEdgeColor; else Color=Obj.Color; end

[data,loc_shape]=XMLaddNode('Shape',data,loc_scene+1);
switch(Obj.Marker(1))
    case {'*','x','+','s','d','o'}
        [data,loc_indexlineset]=XMLaddNode('IndexedLineSet',data,loc_shape+1);
        E(:,3)=0; E=E-1;
        str=array2str(int32(E));
        data=XMLaddProperty('coordIndex',str,data);
        
        data=XMLaddNode('Coordinate',data,loc_indexlineset+1);
        Vr=V;
        Vr(:,1)=Vr(:,1)*data.tags.flip(1);
        Vr(:,2)=Vr(:,2)*data.tags.flip(2);
        Vr(:,3)=Vr(:,3)*data.tags.flip(3);
        str=array2str(Vr);
        data=XMLaddProperty('point',str,data);
    case '.'
        [data,loc_indexfaceset]=XMLaddNode('IndexedFaceSet',data,loc_shape+1);
        data=XMLaddProperty('solid','true',data);
        F(:,4)=0; F=F-1;
        strCoord=array2str(int32(F),',');
        data=XMLaddProperty('coordIndex',strCoord,data);
        data=XMLaddNode('Coordinate',data,loc_indexfaceset+1);
        Vr=V;
        Vr(:,1)=Vr(:,1)*data.tags.flip(1);
        Vr(:,2)=Vr(:,2)*data.tags.flip(2);
        Vr(:,3)=Vr(:,3)*data.tags.flip(3);
        str=array2str(Vr);
        data=XMLaddProperty('point',str,data);
    otherwise
end

[data,loc_appearance]=XMLaddNode('Appearance',data,loc_shape+1);
data=XMLaddNode('Material',data,loc_appearance+1);
if(~isempty(Obj.Color))
    str=sprintf('%4.4f %4.4f %4.4f',Color(1),Color(2),Color(3));
    data=XMLaddProperty('emissiveColor',str,data);
end

function [E,V]=addlinepiece(x1,x2,y1,y2,z1,z2,E,V)
V1=[[x1(:);x2(:)],[y1(:);y2(:)],[z1(:);z2(:)]];
E1=[(1:size(x1,1))' (1:size(x1,1))'+size(x1,1)];
if(nargin<7)
    E=E1; V=V1;
else
    E=[E;E1+size(V,1)];
    V=[V;V1];
end



