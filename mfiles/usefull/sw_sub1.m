function out = sw_sub1(inp, varargin)
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

if nargin == 1
    dnum = 1;
else
    dnum = varargin{1};
end

if isa(inp,'sym')
    
    symVar = symvar(inp);
    
    if ~isempty(symVar)
        out = double(subs(inp,symVar,dnum*ones(1,numel(symVar))));
    else
        out = double(inp);
    end
else
    out = inp;
end

end