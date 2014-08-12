function out = sw_sub1(inp)
% converts symbolic variables into double by substituting 1 for every symbol
%
% out = SW_SUB1(inp)
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
    
symVar = symvar(inp);

if ~isempty(symVar)
    out = double(subs(inp,symVar,ones(1,numel(symVar))));
else
    out = double(inp);
end
else
    out = inp;
end

end