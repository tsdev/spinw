function [symOp, symInfo] = operator(sym, fid)
% generates all symmetry elements from given space group
% 
% ### Syntax
% 
% `[symOp, symInfo] = swsym.operator(sym)`
% 
% `[symOp, symInfo] = swsym.operator(sym,fid)`
%
% ### Description
% 
% `[symOp, symInfo] = swsym.operator(sym)` generates *all* symmetry
% elements from a given set of generators. It also accepts space group
% labels or space group index or string of symmetry operators.
% 
% ### Input Arguments
% 
% `sym`
% : Line index in the [symmetry.dat] file or string of the
%   symmetry operators or matrix of symmetry generators with dimensions of
%   $[3\times 4\times n_{op}]$. For example: `sym = 'P n m a'`.
% 
% `fid`
% : If non-zero, the symmetry operators will be printed to the file
%   identified by `fid`, the following values are valid:
%   * `0`   no printed output (default),
%   * `1`   standard output (Command Line),
%   * `fid` text file opened before using `fid = fopen(path)`.
% 
% ### Output Arguments
% 
% `symOp`
% : All the symmetry elements in a matrix with dimensions of $[3\times
%   4\times n_{op}]$.
%
% `symInfo`
% : Structure that contains additional information about the space 
%   group with the following fields:
%   * `name`    Name of the space group, if the `swsym.generator`
%               function is called with no input, name stores the name of
%               all space groups from [symmetry.dat] file in a cell.
%   * `str`     The string of the symmetry operations.
%   * `num`     The line index in the [symmetry.dat] file.
% 
% ### See Also
% 
% [swsym.generator]
%

if nargin == 0
    help swsym.operator
    return
end

if nargin < 2
    fid = 0;
end

[genOp, symInfo] = swsym.generator(sym);

nGen  = size(genOp,3);
symOp = [eye(3) zeros(3,1)];

% use factor 12 for making integer translation values
genOp(:,4,:) = round(genOp(:,4,:)*12);

% stores the order of the generators
P = zeros(1,nGen);

% generate all symmetry elements
for ii = 1:nGen
    R0 = genOp(:,1:3,ii);
    T0 = genOp(:,4,ii);
    % order of the symmetry operator
    P(ii) = swsym.oporder([R0 T0/12])-1;
    
    R = eye(3);
    T = zeros(3,1);
    
    for jj = 1:P(ii)
        R = R0*R;
        T = R0*T + T0;
        
        nSym = size(symOp,3);
        for kk = 1:nSym
            RS = R*symOp(:,1:3,kk);
            TS = mod(round(R*symOp(:,4,kk)+T),12);
            
            idxR = permute(sumn(abs(bsxfun(@minus,symOp(:,1:3,:),RS)),[1 2]),[3 1 2]);
            idxT = permute(any(bsxfun(@minus,symOp(:,4,:),TS),1),[3 1 2]);
            
            % adds new operator to the list if it differs from all
            if all(idxR | idxT)
                symOp = cat(3,symOp,[RS TS]);
            end
            
        end
    end
end

% convert back translations to fractional numbers
symOp(:,4,:) = symOp(:,4,:)/12;

% print symmetry elements
if fid~=0
    fprintf(fid, 'General coordinates of space group: %s\n',symInfo.name);
    
    sStr = strtrim(strsplit(swsym.str(symOp),';')');
    for ii = 1:numel(sStr)
        fprintf(fid,'(%02d) %s\n',ii,sStr{ii});
    end
end

end