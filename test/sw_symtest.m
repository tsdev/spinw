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
    % check that the number of generators agree
    good(ii) = size(Rg,3) == size(R0,3);
    
end

if ~all(good)
    error('Symmetry generator calculator test failed!')
end














