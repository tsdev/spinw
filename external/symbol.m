function S = symbol(sName)
% returns the character corresponding to the given symbol name
%
% S = SYMBOL(sName)
%
% Input:
%
% sName     String, name of the symbol. For example: 'alpha', 'Angstrom',
%           etc. The input is case sensitive.
%
% Output:
%
% S         A char type variable containing the symbol.
%

nList = {'angstrom' 'copy' 'reg' 'deg' 'pm' 'square' 'cube' 'cross' 'par'...
    'perp' 'int' 'alpha','beta','gamma' 'delta' 'epsilon' 'zeta' 'eta'...
    'theta' 'iota' 'kappa' 'lambda' 'mu' 'nu' 'xi' 'omicron' 'pi' 'rho'...
    'ssigma' 'sigma' 'tau' 'upsilon' 'phi' 'chi' 'psi' 'omega'};

nList = [nList cellfun(@(C)[upper(C(1)) C(2:end)],nList([12:28 30:end]),'UniformOutput',false)];

cList = char([197 169 174 176:179 215 449 10178 8747 945:969 913:929 931:937]);

nL = numel(nList);

if nargin == 0
    help symbol
    fprintf('List of symbols:\n')
    for ii = 1:(nL/2)
        fprintf('%10s%5c%10s%5c\n',nList{ii},cList(ii),nList{ii+nL/2},cList(ii+nL/2));
    end
    return
end

idx = find(cellfun(@(C)~isempty(C)&&(C(1)==1),strfind(nList,sName)),1,'first');

if isempty(idx)
    error('symbol:WrongInput','Symbol with the given name doesn''t exists!')
end

S = cList(idx);

end