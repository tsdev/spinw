function data = addtext(data,loc_scene,Obj)

[data,loc_transform]=XMLaddNode('Transform',data,loc_scene+1);
Pr=Obj.Position.*data.tags.flip;
str=sprintf('%4.4f %4.4f %4.4f',Pr(1),Pr(2),Pr(3));
data=XMLaddProperty('translation',str,data);

fs=Obj.FontSize*(1/data.tags.maxsize);
%if(data.tags.xhtml),  fs=fs/10; else fs=fs/500;  end
if(data.tags.xhtml),  fs=fs/60; else fs=fs/500;  end

str=sprintf('%4.4f %4.4f %4.4f',fs,fs,fs);
data=XMLaddProperty('scale',str,data);

[data,loc_shape]=XMLaddNode('Shape',data,loc_transform+1);

[data,loc_appearance]=XMLaddNode('Appearance',data,loc_shape+1);
[data]=XMLaddNode('Material',data,loc_appearance+1);

str=sprintf('%4.4f %4.4f %4.4f',Obj.Color(1),Obj.Color(2),Obj.Color(3));
data=XMLaddProperty('emissiveColor',str,data);

[data,loc_text]=XMLaddNode('Text',data,loc_shape+1);
data=XMLaddProperty('string',Obj.String,data);
data=XMLaddProperty('solid','false',data);

[data]=XMLaddNode('FontStyle',data,loc_text+1);
data=XMLaddProperty('family','Times',data);
data=XMLaddProperty('size',num2str(16),data);
data=XMLaddProperty('solid','false',data);

end