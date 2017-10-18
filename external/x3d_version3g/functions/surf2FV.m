function [F,V,Cface,Cedge,E,T]=surf2FV(Obj)
    XData=Obj.XData;
    YData=Obj.YData;
    ZData=Obj.ZData;
    if(size(XData,1)==1)
        XData=repmat(XData,size(ZData,1),1);
    end
    if(size(YData,2)==1)
        YData=repmat(YData,1,size(ZData,2));
    end
    [Tx,Ty]=ndgrid(linspace(1,0,size(ZData,1)),linspace(0,1,size(ZData,2)));
    T=[Ty(:) Tx(:)];
    V=[XData(:) YData(:) ZData(:)];
    I=reshape(1:numel(ZData),size(ZData));
    I1=I(1:end-1,1:end-1);
    I2=I(2:end  ,1:end-1);
    I3=I(1:end-1,2:end);
    I4=I(2:end  ,2:end);
    F=[I1(:) I3(:) I4(:); I2(:) I1(:) I4(:)];
    E=[I1(:) I2(:); I2(:) I4(:); I4(:) I3(:);I3(:) I1(:)];
    
    C=Obj.CData;
    if(ischar(Obj.FaceColor)), fc=Obj.FaceColor; else fc=''; end
    switch(fc)
        case 'flat'
            C1=C(1:end-1,1:end-1,:);
            C2=C(2:end  ,1:end-1,:);
            C3=C(1:end-1,2:end,:);
            C4=C(2:end  ,2:end,:);
            C=(C1+C2+C3+C4)/4;
            C=reshape(C,size(F,1)/2,[]);
            C=[C;C];
            C=reshape(C(:),size(F,1),[]);
        case 'interp'
            C=reshape(C(:),size(V,1),[]);
        otherwise
            C=[];
    end
    Cface=C;
    
    C=Obj.CData;
    if(ischar(Obj.EdgeColor)), fc=Obj.EdgeColor; else fc=''; end
    switch(fc)
        case 'flat'
            C1=C(1:end-1,1:end-1,:);
            C2=C(2:end  ,1:end-1,:);
            C3=C(1:end-1,2:end,:);
            C4=C(2:end  ,2:end,:);
            C=(C1+C2+C3+C4)/4;
            C=reshape(C,size(E,1)/4,[]);
            C=[C;C;C;C];
            C=reshape(C(:),size(E,1),[]);
        case 'interp'
            C=reshape(C(:),size(V,1),[]);
        otherwise
            C=[];
    end
    Cedge=C;
    



