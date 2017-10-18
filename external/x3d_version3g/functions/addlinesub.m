function data=addlinesub(data,loc_scene,Obj)
if(~isfield(Obj,'EdgeColor'))
    Obj.EdgeColor=Obj.Color;
end

if(ischar(Obj.EdgeColor)&&strcmpi(Obj.EdgeColor,'none'))
    return;
end

E=Obj.E;
V=Obj.V;

if(~strcmpi(Obj.LineStyle,'none'))
    [data,loc_shape]=XMLaddNode('Shape',data,loc_scene+1);
    [data,loc_indexlineset]=XMLaddNode('IndexedLineSet',data,loc_shape+1);
    str=array2str(int32(E));
    data=XMLaddProperty('coordIndex',str,data);
    
    
    if(isfield(Obj,'EdgeVertexCData')&&(~isempty(Obj.EdgeVertexCData)))
        if(size(Obj.EdgeVertexCData,2)==3)
            cmap = Obj.EdgeVertexCData;
            imap = 0:(size(Obj.EdgeVertexCData,1)-1);
        else
            cmap = colormap;
            imap = double(Obj.EdgeVertexCData(:)'-1);
        end
        
        switch(Obj.CDataMapping)
            case 'scaled'
                imap=imap-min(imap(:));
                imap=imap.*((size(cmap,1)-1)/max(imap(:)));
            case 'direct'
                imap(imap<0)=0;
                imap(imap>(size(cmap,1)-1))=size(cmap,1)-1;
            otherwise
        end
        imap=round(imap);
        
        switch(Obj.EdgeColor)
            case 'flat'
                if(size(imap,2)==size(Obj.Vertices,1))
                    imap=imap(Obj.E(:,1));
                    imap=imap(:)';
                end
                data=XMLaddProperty('colorPerVertex','false',data);
                str=array2str(int32([imap -1]));
                data=XMLaddProperty('colorIndex',str,data);
                
            case 'interp'
                if(size(imap,2)==size(Obj.Vertices,1))
                    data=XMLaddProperty('colorPerVertex','true',data);
                else
                    error('figur2xhtml:err','interp number vertex not equal to color');
                end
            otherwise
        end
    end
    
    data=XMLaddNode('Coordinate',data,loc_indexlineset+1);
    V(:,1)=V(:,1)*data.tags.flip(1);
    V(:,2)=V(:,2)*data.tags.flip(2);
    V(:,3)=V(:,3)*data.tags.flip(3);
    str=array2str(V);
    data=XMLaddProperty('point',str,data);
    
    if(isfield(Obj,'EdgeVertexCData')&&(~isempty(Obj.EdgeVertexCData)))
        data=XMLaddNode('Color',data,loc_indexlineset+1);
        str=array2str(cmap);
        data=XMLaddProperty('color',str,data);
    end
    
    
    [data,loc_appearance]=XMLaddNode('Appearance',data,loc_shape+1);
    data=XMLaddNode('Material',data,loc_appearance+1);
    if(~isempty(Obj.Color))
        str=sprintf('%4.4f %4.4f %4.4f',Obj.Color(1),Obj.Color(2),Obj.Color(3));
        data=XMLaddProperty('emissiveColor',str,data);
    end
    data=XMLaddNode('LineProperties',data,loc_appearance+1);
    data=XMLaddProperty('applied','true',data);
    switch (Obj.LineStyle)
        case '-'
            str='1';
        case ':'
            str='3';
        case '-.'
            str='4';
        case '--'
            str='2';
        otherwise
            str='12';
    end
    data=XMLaddProperty('linetype',str,data);
    str=sprintf('%4.4f',Obj.LineWidth);
    data=XMLaddProperty('linewidthScaleFactor',str,data);
    data=XMLaddProperty('containerField','lineProperties',data);
end

if(~strcmpi(Obj.Marker,'none')),
    data=addlinemarker(data,loc_scene,Obj);
end