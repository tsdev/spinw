function str = sw_rmspace(str)
% remove leading spaces from a cell of strings

% remove empty lines
str(cellfun(@(C)isempty(C),str)) = [];
% lines that begin with space

if all(cellfun(@(C)isempty(C),regexp(str,'^\S','once')))
    sIdx = cellfun(@(C)C(1)==' ' && any(C~=' '),str);
    % minimum space
    nSpace = min(cellfun(@(C)find(diff(C==' '),1,'first'),str(sIdx)));
    % remove leading spaces
    str(sIdx) = cellfun(@(C)C(nSpace+1:end),str(sIdx),'UniformOutput',false);
end

end