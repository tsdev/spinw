function data=addlight(data,loc_scene,Obj)
data=XMLaddNode('PointLight',data,loc_scene+1);
if(~isempty(Obj.Color))
        str=sprintf('%4.4f %4.4f %4.4f',Obj.Color(1),Obj.Color(2),Obj.Color(3));
        data=XMLaddProperty('color',str,data);
end
if(strcmpi(Obj.Style,'local'))
    data=XMLaddProperty('global','false',data);
else
    data=XMLaddProperty('global','true',data);
end
Pr=Obj.Position.*data.tags.flip;
str=sprintf('%4.4f %4.4f %4.4f',Pr(1),Pr(2),Pr(3));
data=XMLaddProperty('location',str,data);
data=XMLaddProperty('intensity','1',data);
data=XMLaddProperty('radius','15',data);
data=XMLaddProperty('ambientIntensity','1',data);
