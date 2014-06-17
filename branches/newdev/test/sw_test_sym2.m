function [Res, errMsg] = sw_test_sym2(tol)

% no error
Res = 0;
errMsg = [];

try
    %% RUN BY HAND
    mnco3 = sw;
    mnco3.genlattice('lat_const',[4.773 4.773 16.642],'angle',[90 90 120]*pi/180,'sym','R -3 c');
    %mnco3.genlattice('lat_const',[4.773 4.773 16.642],'angle',[90 90 120]*pi/180,'sym',167);
    
    mnco3.addatom('label',{'Mn'},'r',[0      0   0]','S',1,'color',[0;0;255]);
    mnco3.addatom('label',{'C'},'r', [0      0 1/4]','S',0,'color',[128;128;128]);
    mnco3.addatom('label',{'O'},'r', [0.2695 0 1/4]','S',0,'color',[255;0;0]);
    
    mnco3.gencoupling('maxdistance',10);
    mnco3.addmatrix('label',{'J1'},'color',[255;0;0]);
    mnco3.addmatrix('label',{'J2'},'color',[0;255;0]);
    mnco3.addmatrix('label',{'J3'},'color',[0;0;255]);
    
    mnco3.addcoupling('J1',1);
    mnco3.addcoupling('J2',1);
    %mnco3.addcoupling('J3',3);
    
    mnco3.setmatrix('label','J1','pref',{[1 0 0]});
    mnco3.setmatrix('label','J2','pref',{[0 1 0]});

    hFig = plot(mnco3,'range',[2 2 1],'pNonMagAtom',false);
    %%
    close(hFig);
    
catch errMsg
    % code throws error
    Res = 1;
    return;
end

end