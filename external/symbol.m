function S = symbol(sName,noError)
% returns the character corresponding to the given symbol name
%
% S = SYMBOL(sName)
%
% Input:
%
% sName     String, name of the symbol. For example: 'alpha', 'Angstrom',
%           etc. The input is case sensitive. Alternatively it can be a
%           text that contains name of symbols with \\ before them, for
%           example: "The length of the bond is 3.2 \\a.".
%
% Output:
%
% S         A char type variable containing the symbol.
%

if any(sName == '\')
    % exchange symbols referenced by '\\symbolname'
    S = regexprep(sName,'\\\\([\w\^-]+)','${symbol($1,2)}');
    return
end

nList = {'hbar' 'angstrom' 'copy' 'reg' 'deg' 'pm' 'square' 'cube' 'cross' 'par'...
    'perp' 'int' 'Euro' 'alpha','beta','gamma' 'delta' 'epsilon' 'zeta' 'eta'...
    'theta' 'iota' 'kappa' 'lambda' 'mu' 'nu' 'xi' 'omicron' 'pi' 'rho'...
    'ssigma' 'sigma' 'tau' 'upsilon' 'phi' 'chi' 'psi' 'omega'};

nList = [nList cellfun(@(C)[upper(C(1)) C(2:end)],nList([14:28 30:end]),'UniformOutput',false)];


cList = [295 197 169 174 176:179 215 449 10178 8747 8364 945:969 913:929 931:937];

% add other symbols
nList = [nList {'skull' 'sun' 'moon' 'ok'   'Angstrom' '^-' '_-' 'bra' 'ket'}];
cList = [cList  9760    9788  9789   10004  197        8315 8331 10216 10217 ];

nList = [nList {'leq' 'geq' 'equiv' 'll' 'gg' 'propto'}];
cList = [cList  8818  8819  8801    8810 8811  8733    ];



nList = [nList arrayfun(@(E)['^' E],num2str((0:9)'),'UniformOutput',false)' arrayfun(@(E)['_' E],num2str((0:9)'),'UniformOutput',false)'];
cList = [cList 8304 185 178:179 8308:8313 8320:8329];

% convert to characters
cList = char(cList);
nL    = numel(nList);

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
    if nargin>1 && any(strcmp(noError,{'copy' 'c'}))
        clipboard('copy',S);
    end
elseif isempty(idx) && nargin>1 && noError==2
    S = sName;
elseif isempty(idx) && nargin>1 && noError
    S = [':::' sName ':::'];
    warning(['The following symbol does not exists: ''' sName '''\n'])
else
    error('symbol:WrongInput','Symbol ''%s'' with the given name doesn''t exists!',sName)
end

end