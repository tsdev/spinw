function r = sw_cmod(r, tol)
% modulo one with tolerance
%
% r = SW_CMOD(r, tol)
%
% It calculates modulo one with tolerance, numbers larger than 1-epsilon >
% 1-tol will be converted to -epsilon.
%
% See also MOD.
%

if nargin == 0
    help sw_cmod
    return
end

r = mod(r,1);

r(r > 1-tol) = r(r > 1-tol)-1;

end