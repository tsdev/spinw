function B  = changem(A,newval,oldval)
% vectorized version of changem
%
% B = CHANGEM(A,newval,oldval)
%
% The function changes the old values in A to the new values pair wise.
%

B = A;
[valid,id] = max(bsxfun(@eq,A(:),oldval(:).'),[],2);
B(valid)   = newval(id(valid));

% for ii = 1:numel(newval)
%     A(A==oldval(ii)) = newval(ii);
% end

end