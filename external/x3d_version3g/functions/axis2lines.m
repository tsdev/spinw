function axis2lines()
view(3);
axis equal; drawnow;
Obj=get(gca);
axis off;
hold on;
r=0.04*max(max(Obj.XLim(2)-Obj.XLim(1),Obj.YLim(2)-Obj.YLim(1)),Obj.ZLim(2)-Obj.ZLim(1));
plot3([Obj.XLim(1) Obj.XLim(2)],[Obj.YLim(1) Obj.YLim(1)],[Obj.ZLim(1) Obj.ZLim(1)],'k');
x=Obj.XTick;
y=ones(size(Obj.XTick))*Obj.YLim(1);
z=ones(size(Obj.XTick))*Obj.ZLim(1);
plot3([x;x],[y-r;y],[z;z],'k');
for i=1:length(Obj.XTick)
    text(x(i),y(i)-2*r,z(i),num2str(Obj.XTick(i)),'FontSize',8);
end

plot3([Obj.XLim(1) Obj.XLim(1)],[Obj.YLim(1) Obj.YLim(2)],[Obj.ZLim(1) Obj.ZLim(1)],'k');
x=ones(size(Obj.YTick))*Obj.XLim(1);
y=Obj.YTick;
z=ones(size(Obj.YTick))*Obj.ZLim(1);
plot3([x-r;x],[y;y],[z;z],'k');
for i=1:length(Obj.YTick)
    text(x(i)-2*r,y(i),z(i),num2str(Obj.YTick(i)),'FontSize',8);
end

plot3([Obj.XLim(1) Obj.XLim(1)],[Obj.YLim(2) Obj.YLim(2)],[Obj.ZLim(1) Obj.ZLim(2)],'k');
x=ones(size(Obj.ZTick))*Obj.XLim(1);
y=ones(size(Obj.ZTick))*Obj.YLim(2);
z=Obj.ZTick;
plot3([x-r;x],[y;y],[z;z],'k');
for i=1:length(Obj.ZTick)
    text(x(i)-2*r,y(i),z(i),num2str(Obj.ZTick(i)),'FontSize',8);
end
view(3);
 
FV.vertices=[0 0 0; 0 0 1; 0 1 0; 0 1 1; 1 0 0; 1 0 1; 1 1 0; 1 1 1];
FV.faces=[1 2 4; 1 4 3;5 8 6; 5 7 8; 1 5 2; 5 6 2; 3 4 7; 7 4 8; 1 3 5; 5 3 7; 2 6 4; 6 8 4];
FV.faces=FV.faces(:,[3 2 1]);
FV.vertices(:,1)=FV.vertices(:,1)*(Obj.XLim(2)-Obj.XLim(1))+Obj.XLim(1);
FV.vertices(:,2)=FV.vertices(:,2)*(Obj.YLim(2)-Obj.YLim(1))+Obj.YLim(1);
FV.vertices(:,3)=FV.vertices(:,3)*(Obj.ZLim(2)-Obj.ZLim(1))+Obj.ZLim(1);
patch(FV,'facecolor',[1 1 1],'facealpha',0.8,'edgecolor','none');
