function g = gcdv(v)
% calculates the greates common divisor of a set of numbers
%
% G = VECGCD(V)
%
% is the greatest common divisor of the elements of the integer vector V. V
% must be a row or a column vector of integers. We define gcd([]) = 1 and
% gcd([0 0 ... 0]) = 0.
%

[r,c] = size(v);
if (r>1) && (c>1)
    error('v must be a row or a column vector');
end

% vecgcd([]) = 1
if isempty(v)
    g = 1;
    return
end

% vecgcd([0 0 ... 0]) = 0
if all(v==0)
    g = 0;
    return
end


% I want a to work with a row vector
if c==1
    v = v';
end

% remove the zeros and sort
v = v(v>0);
v = sort(abs(v),2,'descend');

i = 1;
j = length(v);
while i < j
    x = mod(v(i), v(j));
    x = min([x, v(j) - x]);
    % x|v(i), x|v(j), g|x
    % and x is smaller than min[v]
    
    i = i + 1; % "remove" v(i)
    if x > 0 % if x > 0 append it to v
        j = j + 1;
        v = [v, x]; %#ok<AGROW>
    end
end

% x = v(i)
% x divides v(1), v(2) , ..., v(end)
% g|x
% hence g = x.

g = v(i);

end