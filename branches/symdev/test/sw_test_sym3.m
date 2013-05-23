function [Res, errMsg] = sw_test_sym3(tol)

% no error
Res = 0;
errMsg = [];

try
    %% RUN BY HAND
    
    % LuVO3 intermediate temperature Pbnm space group
    
    a = 5.2821;
    b = 5.6144;
    c = 7.5283;
       
    sw_addsym('x+1/2,-y+1/2,-z; -x,-y,z+1/2; -x,-y,-z','P b n m');
    
    luvo = sw;
    luvo.genlattice('lat_const',[a b c],'sym','P b n m');
    
    % V1 atom
    luvo.addatom('r',[1/2 0 0],'label','V','S',1,'color',[128; 128; 128]);
    % O1 atom
    luvo.addatom('r',[0.11560;0.44821;1/4],'label','O1','S',0,'color',[255; 0; 0]);
    % O2 atom
    luvo.addatom('r',[0.68170;0.30500;0.06540],'label','O2','S',0,'color',[255; 0; 0]);
    
    luvo.gencoupling('maxdistance',8);
    
    luvo.addmatrix('label',{'J1a'});
    luvo.addmatrix('label',{'J1b'});
    luvo.addmatrix('label',{'Jab'},'color',[0;0;255]);
    
    luvo.addcoupling('J1a',1);
    luvo.addcoupling('J1b',1);
    
    luvo.setmatrix('label','J1a','pref',{[1 0 0]});
    luvo.setmatrix('label','J1b','pref',{[0 1 0]});
    
    hFig = plot(luvo,'pNonMagAtom',false,'pZeroCoupling',true);
    
    %%
    
    close(hFig);

catch errMsg
    % code throw error
    Res = 1;
    return;
end

try
    
    cMat = luvo.getmatrix('coupling_idx',1);
    sMat = (cMat-permute(cMat,[2 1 3])).^2;
    sMat = permute(sum(sum(sMat)),[2 3 1]) > tol;
    
    if sum(sMat) ~= 2
        error('sw_test_sym1:WrongMat','Wrong number of allowed antisymmetric coupling matrix!');
    end
    
    % test the generated DM interactions
    gMat = reshape(luvo.intmatrix.all(end-8:end,:),3,3,[]);
    % (- - + + - + + -)
    if any(abs(permute(cat(3,gMat(2,3,1:4),gMat(1,3,5:8)),[2,3,1]) - [-1 -1 1 1 -1 1 1 -1])>tol)
        error('sw_test_sym1:WrongMat','Wrong values of the generated antisymmetric coupling matrix!');
    end
    if any(abs(permute(cat(3,gMat(3,2,1:4),gMat(3,1,5:8)),[2,3,1]) - [1 1 -1 -1 1 -1 -1 1])>tol)
        error('sw_test_sym1:WrongMat','Wrong values of the generated antisymmetric coupling matrix!');
    end
    
    
catch errMsg
    % wrong result
    Res = 2;
end

end