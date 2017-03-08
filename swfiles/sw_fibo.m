function [F,F1] = sw_fibo(Fmax)
% returns the last two Fibonacci number smaller or equal to the
% given number
%
% [Flast Fprev] = sw_fibo(Fmax)
%

num = [0 0 1];

while num(end)<Fmax
    num(end+1) = sum(num(end+[-1 0])); %#ok<AGROW>
end

if num(end) == Fmax
    F = num(end);
    F1 = num(end-1);
else
    F = num(end-1);
    F1 = num(end-2);
end
end