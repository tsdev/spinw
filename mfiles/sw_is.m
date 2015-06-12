function out = sw_is(inp)
% check if symbolic logical expression is always valid
%
% WARNING WARNING WARNING WARNING WARNING WARNING WARNING WARNING WARNING 
% 
% This is not a strict test, it only checks whether the expression contains
% a symbolic variable, then the output is false. If the symbolic expression
% contains only constants, it returns the right logical value.
% For example:
% sw_is(sym('1')==1)
% ans = true
%
% sw_is(sym('a','positive')>-1)
% ans = false
%
% WARNING WARNING WARNING WARNING WARNING WARNING WARNING WARNING WARNING 
%
% out = SW_IS(inp)
%
% Input:
%
% inp   Any symbolic/double type matrix.
%
% Output:
%
% out   Double type output with the same dimensions as the input.
%

if isa(inp,'sym')
    % find elements without symbolic variable
    selNum = arrayfun(@(x)isempty(symvar(x)),inp);
    
    out = false(size(inp));
    if any(selNum)
        out(selNum) = isAlways(inp(selNum));
    end

else
    out = logical(inp);
end

end