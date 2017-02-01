function A  = changem(A,newval,oldval)
% vectorized version of changem
%
% B = CHANGEM(A,newval,oldval)
%
% The function changes the old values in A to the new values pair wise.
%

[idx1, idx2] = ismember(A,oldval);
A(idx1) = newval(idx2(idx1));

end