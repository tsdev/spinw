function S = symbol(sName,noError)
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

nList = {'hbar' 'angstrom' 'copy' 'reg' 'deg' 'pm' 'square' 'cube' 'cross' 'par'...
    'perp' 'int' 'Euro' 'alpha','beta','gamma' 'delta' 'epsilon' 'zeta' 'eta'...
    'theta' 'iota' 'kappa' 'lambda' 'mu' 'nu' 'xi' 'omicron' 'pi' 'rho'...
    'ssigma' 'sigma' 'tau' 'upsilon' 'phi' 'chi' 'psi' 'omega'};

nList = [nList cellfun(@(C)[upper(C(1)) C(2:end)],nList([14:28 30:end]),'UniformOutput',false)];


cList = char([295 197 169 174 176:179 215 449 10178 8747 8364 945:969 913:929 931:937]);

% add other symbols
nList = [nList      {'skull' 'sun' 'moon' 'ok'   'Angstrom'}];
cList = [cList char([9760    9788  9789   10004  197      ])];

nL = numel(nList);

if nargin == 0
    help symbol
    fprintf('List of symbols:\n')
    for ii = 1:floor(nL/2)
        fprintf('%10s%5c%10s%5c\n',nList{ii},cList(ii),nList{ii+floor(nL/2)},cList(ii+floor(nL/2)));
    end
    if ceil(nL/2)~=floor(nL/2)
        fprintf('%10s%5c\n',nList{end},cList(end));
    end
    return
end

idx = find(cellfun(@(C)~isempty(C)&&(C(1)==1),strfind(nList,sName)),1,'first');

if ~isempty(idx)
    S = cList(idx);
elseif isempty(idx) && nargin>1 && noError
    %S = [':::' sName ':::'];
    S = sName;
    fprintf(['The following symbol does not exists: ''' sName '''\n'])
else
    error('symbol:WrongInput','Symbol ''%s'' with the given name doesn''t exists!',sName)
end

end