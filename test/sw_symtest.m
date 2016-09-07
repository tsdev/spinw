%% test symmetry string generator
good = false(1,230);

for ii = 1:230
    [R,T,~,strSym0] = sw_gensym(ii);
    strSym0 = strtrim(strSym0);
    
    strSym = sw_gensymstr(R,T);
    
    % check equality
    if numel(strSym0) == numel(strSym)
        good(ii) = all(strSym0==strSym);
    end
    
end

if ~all(good)
    error('Symmetry string generator test failed!')
end

%% test generator function

good = false(1,230);

for ii = 1:230
    [Rg,Tg] = sw_gensym(ii);
    [R,T]   = sw_gencoord(ii);
    [R0, T0] = sw_symgetgen(R,T);
    % check that the number of generators are smaller or equal 
    good(ii) = size(Rg,3) >= size(R0,3);
    
end

[R1,T1]   = sw_gencoord({R0 T0});
[R2,T2]   = sw_gencoord({Rg(:,:,[1 3 4]) Tg(:,[1 3 4])});

if ~all(good)
    error('Symmetry generator calculator test failed!')
end

% %%
% 
% Rs =[reshape(R,9,[]); T];
% R1s=[reshape(R1,9,[]);T1];















