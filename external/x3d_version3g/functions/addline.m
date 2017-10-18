function data=addline(data,loc_scene,Obj)
    if(length(Obj.XData)==length(Obj.ZData))
        V=[Obj.XData(:) Obj.YData(:) Obj.ZData(:)];
    else
        V=[Obj.XData(:) Obj.YData(:)]; V(:,3)=0;
    end
    
    Obj.E=0:size(V,1)-1;
    Obj.V=V;
    data=addlinesub(data,loc_scene,Obj);
   
