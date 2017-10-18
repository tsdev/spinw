function data=addmesh(data,loc_scene,Obj)
data.tags.numobjects=data.tags.numobjects+1;

[data,loc_shape]=XMLaddNode('Shape',data,loc_scene+1);
[data,loc_appearance]=XMLaddNode('Appearance',data,loc_shape+1);
data=XMLaddNode('Material',data,loc_appearance+1);

if(isempty(Obj.FaceColor)), Obj.FaceColor=[1 0 0]; end

if(data.tags.xhtml)
    data=XMLaddProperty('id',['mshape' num2str(data.tags.numobjects-1)],data);
end

if((~ischar(Obj.FaceAlpha))&&(~isempty(Obj.FaceAlpha)))
    str=sprintf('%4.4f ',1-Obj.FaceAlpha);
    data=XMLaddProperty('transparency',str,data);
end

%str=sprintf('%4.4f %4.4f %4.4f',1,1,1);
%data=XMLaddProperty('specularColor',str,data);
%data=XMLaddProperty('ambientIntensity','1',data);
%data=XMLaddProperty('shininess','1',data);

if(isnumeric(Obj.FaceColor))
    str=sprintf('%4.4f %4.4f %4.4f',Obj.FaceColor(1),Obj.FaceColor(2),Obj.FaceColor(3));
    data=XMLaddProperty('diffuseColor',str,data);
else
    switch(Obj.FaceColor)
        case 'texturemap'
            [data,loc_texture]=XMLaddNode('Texture',data,loc_appearance+1);
            data=XMLaddProperty('hideChildren','false',data);
            filename=[data.tags.filename 't' num2str(data.tags.numobjects) '.png'];
            I=Obj.CData;
            if(size(I,3)==1)
                cmap = colormap;
                switch(Obj.CDataMapping)
                    case 'scaled'
                        I=I-min(I(:));
                        I=round(I.*((size(cmap,1)-1)/max(I(:))));
                        Ic=cmap(I(:)+1,:);
                        I=reshape(Ic,[size(I) 3]);
                    case 'direct'
                        I(I<0)=0;
                        I(I>(size(cmap,1)-1))=size(cmap,1)-1;
                        Ic=cmap(I(:)+1,:);
                        I=reshape(Ic,[size(I) 3]);
                    otherwise
                end
            end
            if(ischar(Obj.FaceAlpha)&&strcmpi(Obj.FaceAlpha,'texturemap'))
                A=Obj.AlphaData;
                amap = colormap;
                switch(Obj.CDataMapping)
                    case 'scaled'
                        A=A-min(A(:));
                        A=round(A.*((size(amap,1)-1)/max(A(:))));
                        Ac=amap(A(:)+1,:);
                        A=reshape(Ac,size(A));
                    case 'direct'
                        A(A<0)=0;
                        A(A>(size(amap,1)-1))=size(amap,1)-1;
                        Ac=amap(A(:)+1);
                        A=reshape(Ac,size(A));
                    otherwise
                end
            else
                A=[];
            end
            if(data.tags.xhtml)
				if(data.tags.options.embedimages)
					filename='temptex.png';
				end
                s=2^ceil(log(max(size(I)))/log(2));
                I=imresize(I,[s s]);
                if(~isempty(A)), A=imresize(A,[s s]); end
                if(isempty(A))
                    imwrite(I,[data.tags.folder filename]);
                else
                    imwrite(I,[data.tags.folder filename],'Alpha',A);
                end
                
				if(data.tags.options.embedimages)
					imagebase64=base64_encode([data.tags.folder filename]);
					delete([data.tags.folder filename]);
					data=XMLaddNode('Texture',data,loc_texture+1);
					data=XMLaddProperty('src',['data:image/png;base64,' imagebase64],data);
				else
				    data=XMLaddNode('Texture',data,loc_texture+1);
					data=XMLaddProperty('src',filename,data);
				end
            else
                if(isempty(A))
                    imwrite(I,[data.tags.folder filename]);
                else
                    imwrite(I,[data.tags.folder filename],'Alpha',A);
                end
                data=XMLaddNode('ImageTexture',data,loc_texture+1);
                data=XMLaddProperty('url',filename,data);
            end
            
    end
end


[data,loc_indexfaceset]=XMLaddNode('IndexedFaceSet',data,loc_shape+1);
E1=[Obj.Faces(:,[1 2]);Obj.Faces(:,[2 3]);Obj.Faces(:,[3 1]);Obj.Faces(:,[2 1]);Obj.Faces(:,[3 2]);Obj.Faces(:,[1 3])];
E2=unique(E1,'rows');
solid=size(E1,1)==(size(E2,1)*2);
if(solid)
    data=XMLaddProperty('solid','true',data);
else
    data=XMLaddProperty('solid','false',data);
end
if(data.tags.xhtml&&data.tags.options.interactive)
    data=XMLaddProperty('onclick',['clickalert(' num2str(data.tags.numobjects-1) ');'],data);
end

F=Obj.Faces-1;
F(:,4)=-1;
strCoord=array2str(int32(F),',');

if(~isempty(Obj.FaceVertexCData))
    if(size(Obj.FaceVertexCData,2)==3)
        cmap = Obj.FaceVertexCData;
        imap = 0:(size(Obj.FaceVertexCData,1)-1);
    else
        cmap = colormap;
        imap = double(Obj.FaceVertexCData(:)'-1);
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
    
    switch(Obj.FaceColor)
        case 'flat'
            if(size(imap,2)==size(Obj.Vertices,1))
                imap=imap(Obj.Faces(:,1));
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
            
        case 'texturemap'
        otherwise
    end
end

if(ischar(Obj.FaceColor)&&strcmpi(Obj.FaceColor,'texturemap'))
    data=XMLaddProperty('texCoordIndex',strCoord,data);
end

data=XMLaddProperty('coordIndex',strCoord,data);
data=XMLaddNode('Coordinate',data,loc_indexfaceset+1);
V=Obj.Vertices;
Vr=V;
Vr(:,1)=Vr(:,1)*data.tags.flip(1);
Vr(:,2)=Vr(:,2)*data.tags.flip(2);
Vr(:,3)=Vr(:,3)*data.tags.flip(3);
str=array2str(Vr);
data=XMLaddProperty('point',str,data);

if(ischar(Obj.FaceColor)&&strcmpi(Obj.FaceColor,'texturemap'))
    data=XMLaddNode('TextureCoordinate ',data,loc_indexfaceset+1);
    VT=Obj.TextureVertices;
    str=array2str(VT);
    data=XMLaddProperty('point',str,data);
end



if(~isempty(Obj.FaceVertexCData))
    data=XMLaddNode('Color',data,loc_indexfaceset+1);
    str=array2str(cmap);
    data=XMLaddProperty('color',str,data);
end

if(~isfield(Obj,'E'))
    F(:,4)=F(:,1); F(:,5)=-1;
    Obj.E=F;
else
    Obj.E=Obj.E-1;
    Obj.E(:,3)=-1;
end

Obj.V=V;
Obj.Color=Obj.EdgeColor;
data=addlinesub(data,loc_scene,Obj);
% Lines have the same coordinates as the faces, thus GPU will
% show some line-pieces because they are in front of the face
% duo to some round-off error, and hide some line pieces because
% they are behind the face due to round-off error.
%
% Therefore we add extralines, behind and before the faces
if(isfield(Obj,'VertexNormals'))
    d=mean(max(V,[],1)-min(V,[],1));
    N=reshape(Obj.VertexNormals,size(V));
    L=sqrt(N(:,1).^2+N(:,2).^2+N(:,3).^2)+eps;
    N(:,1)=N(:,1)./L; N(:,2)=N(:,2)./L; N(:,3)=N(:,3)./L;
    Obj.V=V+d*N*0.001;
    data=addlinesub(data,loc_scene,Obj);
    Obj.V=V-d*N*0.001;
    data=addlinesub(data,loc_scene,Obj);
end

