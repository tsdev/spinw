function [Res, errMsg] = sw_test_sym1(tol)
% test symmetry on different simple cyclic space groups

% no error
Res = 0;
errMsg = [];

try
    %% TEST RUN BY HAND
    tetra = sw;
    tetra.genlattice('lat_const',[6 6 5],'sym','P 4','angle',[90 90 90]*pi/180);
    
    %tetra.addatom('r',[1/4 1/4 0]);
    tetra.addatom('r',[0/2 0/2 0]);
    
    tetra.addmatrix('label',{'J1'});
    tetra.addmatrix('label',{'J2'},'color',[0; 255; 0]);
    tetra.addmatrix('label',{'A'},'mat',diag([1 0 0]));
    
    tetra.gencoupling('maxDistance',10);
    tetra.addcoupling('J1',1);
    tetra.addcoupling('J2',2);
    
    
    tetra.addaniso('A',1);
    %tetra.matrix.mat(:,:,3) = [0 1 0;1 0 0;0 0 0];
    tetra.setmatrix('label','A','fid',0,'pref',{[1 1]})
    
    hFig = plot(tetra,'range',[2 2 1]);
    
    %%
    close(hFig);
    
catch errMsg
    
    % code throws error
    Res = 1;
    return;
end

% test whether the results are right
try
    aMat = tetra.getmatrix('aniso_idx',1);
    if size(aMat,3) ~= 2
        error('sw_test_sym1:WrongMat','Wrong number of allowed anisotropy matrix!');
    end
    rightMat = cat(3,diag([0 0 1]),diag([1 1 0]));
    errMat1 = (aMat - rightMat).^2;
    errMat2 = (aMat - rightMat(:,:,[2 1])).^2;
    if (sum(errMat1(:))>tol) && (sum(errMat2(:)) > 1e-5)
        error('sw_test_sym1:WrongMat','Wrong values in the allowed anisotropy matrix!');
    end
    
    cMat = tetra.getmatrix('coupling_idx',1);
    if size(cMat,3) ~= 3
        error('sw_test_sym1:WrongMat','Wrong number of allowed coupling matrix!');
    end
    rightMat = cat(3,diag([0 0 1]),diag([1 1 0]),[0 -1 0;1 0 0;0 0 0]);
    for ii = 1:3
        right = false;
        for jj = 1:3
            errMat = (cMat(:,:,ii) - rightMat(:,:,jj)).^2;
            if sum(errMat(:)) < tol
                right = true;
            end
        end
        if ~right
            error('sw_test_sym1:WrongMat','Wrong values in the allowed coupling matrix!');
        end
    end
    
catch errMsg
    
    % wrong result
    Res = 2;
end

end