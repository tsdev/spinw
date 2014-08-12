function out = sw_always(inp)
% converts symbolic logical expressions into logical expressions
%
% out = SW_ALWAYS(inp)
%
% Input:
%
% inp   Any symbolic/logical type matrix.
%
% Output:
%
% out   Logical output with the same dimensions as the input.
%

if isa(inp,'sym')
    out = isAlways(inp);
else
    out = logical(inp);
end

end