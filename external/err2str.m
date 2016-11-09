function sOut = err2str(num0,err)
% converts value and standard deviation into a string
%
% sOut = ERR2STR(num,{err})
%
% The result is a number with error in brackets in the end. All standard
% deviation is given with the 2 leading digits if the 2 leading digits of
% the standard deviation is smaller than 20. Otherwise the first leading
% digit of the s.d. is given. Function also removes trailing zero to
% improve the quality of the string.
%
% Input:
%
% num   Single number of a 2 element vector, where the second element is
%       the standard deviation.
% err   Standard deviation. If given, this will be used instead of num(2).
%       Optional.
%
% WARNING it gives sometimes wrong results, e.g. err2str([100.1 1])
%

num = num0;

if nargin<2 && numel(num)>1
    err = num(2);
    num = num(1);
elseif nargin<2
    err = 0;
end

err = abs(err);
sgn = sign(num);
num = abs(num);

% get the first digit of the error
eDig = floor(log10(err))-1;
e10 = err/10^eDig;

if e10>=20
    e10 = e10/10;
end

% digits of error
if err>0
    eExp = round(log10(err/e10));
else
    eExp = 0;
end

if err>=10
    if e10>=20
        e10 = round(err,-eDig+1);
    else
        e10 = round(err,-eDig+2);
    end
end

% round the value for the given digits
numR = round(num,-eExp);

% number of digits integer
nDig = floor(log10(numR))+1;

if eExp<0
    if e10 == 10 && nDig>1 && eExp == -1
        nDec = 0;
        e10 = 1;
    else
        nDec = -eExp;
    end
else
    nDec = 0;
end

if nDig<1
    nDig = 1;
end

nInt = nDig+1+nDec;

e10 = round(e10);

s1 = num2str(numR,['%' num2str(nInt) '.' num2str(nDec) 'f']);
s2 = num2str(e10,'%d');

if nDec>0 && s1(end) == '0' && s2(end) == '0'
    s1 = s1(1:(end-1));
    s2 = s2(1:(end-1));
end

if s1(end) == '.'
    s1 = s1(1:(end-1));
end

if err == 0
    sOut = num2str(num0(1));
else
    if sgn > 0
        sOut = [s1 '(' s2 ')'];
    else
        sOut = ['-' s1 '(' s2 ')'];
    end
end

end