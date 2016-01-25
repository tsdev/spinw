function S = symbol(sName)
% returns the character corresponding to the given symbol name
%
% S = SYMBOL(sName)
%
% Input:
%
% sName     String, name of the symbol. For example: 'alpha', 'Angstrom',
%           etc.
%
% Output:
%
% S         A char type variable containing the symbol.
%

nList = {'angstrom','alpha','beta','gamma' 'delta'};

cList = char([197 945:948]);


idx = find(cellfun(@(C)~isempty(C)&&(C(1)==1),strfind(nList,sName)),1,'first');

if isempty(idx)
    error('symbol:WrongInput','Symbol with the given name doesn''t exists!')
end

S = cList(idx);

end