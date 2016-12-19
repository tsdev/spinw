function symStr = str(symOp)
% generates a string equivalent of symmetry operators
%
% symStr = SWSYM.STR(symOp)
%
% Input:
%
% symOp     Symmetry operator with rotations matrices symOp(:,1:3,:) and
%           translation vectors in symOp(:,4,:).
%
% Output:
%
% strSym    String, contains the symmetry operations.
%
% See also SWSYM.ADD, SWSYM.GENERATOR.
%

if nargin == 0
    help swsym.str
    return
end

% number of perators
nOp = size(symOp,3);

% string for negative sign, & is used as a placeholder
sgnStr =  '-+';
% string for new operator
newStr = '&&;';
% string to separate x,y,z component
comStr = ',,&';
% space before each operator
spStr  = '&& ';
% division sign for translation
perStr = repmat('/',[1 nOp*3]);
% each column of R gives the new x1', y1', z1', x2', ...
R = reshape(permute(symOp(:,1:3,:),[2 1 3]),3,[]);
% each column gives the translation dx1, dy1, dz1, dx2, ...
T = symOp(:,4,:);
T = T(:)';
% create the signs for x,y,z
sgnR = sgnStr((R>=0)+1);
% create x,y,z symbols
strR = bsxfun(@times,('xyz')',abs(R));
sgnR(strR==0) = '&';

% remove + signs if there is no previous symbol
sgnR(1,sgnR(1,:)=='+') = '&';
sgnR(2,strR(1,:)==0 & sgnR(2,:)=='+') = '&';
sgnR(3,strR(1,:)==0 & strR(2,:)==0 & sgnR(3,:)=='+') = '&';
strR(strR==0) = '&';
% put all together
strR = [sgnR;strR];
strR = strR([1 4 2 5 3 6],:);
% generate the translations from the nominator, denominator
[Tn,Td] = rat(T);
Td(Tn==0) = 0;
perStr(Tn==0) = '&';
% get the sign of T
sgnT = abs((Tn>=0)+1/2)+1/2;
strS = sgnStr(sgnT);
strS(Tn==0) = '&';
strT = [strS;num2str(Tn')';perStr;num2str(Td')';...
    comStr(repmat(1:3,1,nOp));newStr(repmat(1:3,1,nOp));spStr(repmat(1:3,1,nOp))];
% remove zeros
strT(strT=='0') = '&';
% operator
symStr = [strR;strT];
symStr(symStr=='&') = [];
symStr = symStr(1:(end-2));

end