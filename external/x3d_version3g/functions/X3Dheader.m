function [data,loc_scene]=X3Dheader(data,loc_body,haxis)
	if(data.tags.xhtml)
        [data,loc_x3d]=XMLaddNode('X3D',data,loc_body+1);
		data=XMLaddProperty('xmlns','http://www.web3d.org/specifications/x3d-namespace',data);
		data=XMLaddProperty('showStat','false',data);
		data=XMLaddProperty('showLog','false',data);
		data=XMLaddProperty('x','0px',data);
		data=XMLaddProperty('y','0px',data);
		data=XMLaddProperty('width',[num2str(data.tags.options.width) 'px'],data);
		data=XMLaddProperty('height',[num2str(data.tags.options.height) 'px'],data);
    else
        tags=data.tags;
	    [data,loc_x3d]=XMLaddNode('X3D');
        data.tags=tags;
		data=XMLaddProperty('version','3.1',data);
		data=XMLaddProperty('profile','Immersive',data);
		data=XMLaddProperty('xmlns:xsd','http://www.w3.org/2001/XMLSchema-instance',data);
		data=XMLaddProperty('xsd:noNamespaceSchemaLocation','http://www.web3d.org/specifications/x3d-3.1.xsd',data);
	end
	
	[data,loc_scene]=XMLaddNode('Scene',data,loc_x3d+1);
    data=XMLaddNode('Viewpoint',data,loc_scene+1);
    a=axis(haxis);
    ca=(a(1:2:end)+a(2:2:end))/2;
    da=max(abs(a(1:2:end)-a(2:2:end)));
    if(length(ca)<3), ca(3)=mean(ca); end
    data.tags.maxsize=da;
    flip =[0 0 0];
    if(strcmpi(get(haxis,'XDir'),'normal')), flip(1)=1; else flip(1)=-1; end
    if(strcmpi(get(haxis,'YDir'),'normal')), flip(2)=1; else flip(2)=-1; end
    if(strcmpi(get(haxis,'ZDir'),'normal')), flip(3)=1; else flip(3)=-1; end
    data.tags.flip=flip;
    
    str=sprintf('%4.4f %4.4f %4.4f',flip(1)*ca(1),flip(2)*ca(2),flip(3)*(ca(3)+da));
    data=XMLaddProperty('position',str,data);
    str=sprintf('%4.4f %4.4f %4.4f %4.4f',0,0,1,0);
    data=XMLaddProperty('orientation',str,data);

    data=XMLaddNode('NavigationInfo',data,loc_scene+1);
    if(data.tags.options.headlight)
        data=XMLaddProperty('headlight','true',data);
    else
        data=XMLaddProperty('headlight','false',data);
    end
