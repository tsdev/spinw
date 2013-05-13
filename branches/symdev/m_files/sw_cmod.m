function r = sw_cmod(r, tol)
% modulo one with tolerance, numbers larger than 1-epsilon > 1-tol will be
% converted to -epsilon.
%

r = mod(r,1);

r(r > 1-tol) = r(r > 1-tol)-1;

end
