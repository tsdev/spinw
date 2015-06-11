function V2 = interm(V)
% calculate intermediate point vector
%
% V2 = interm(V)
%
% V2 will have one element less.

if nargin == 0
    help interm
    return
end

V2 = (V(2:end)+V(1:end-1))/2;

end