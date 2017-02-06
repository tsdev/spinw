function [symOp, symInfo] = operator(sym, fid)
% calculates all symmetry operators or general positions for a space group
%
% [symOp, symInfo] = SWSYM.OPERATOR(sym, fid)
%
% Input:
%
% sym           Line index in the symmetry.dat file or string of the
%               symmetry operators or matrix of symmetry generators with
%               dimensions of [3 4 nOp].
%               For example:
%                   sym = 'P b n m';
% fid           For printing the symmetry operators:
%                   0   no printed output (Default)
%                   1   standard output (Command Line)
%                   fid text file opened before with the fid = fopen(path)
%
% Output:
%
% symOp         The rotational part of the symmetry operators, dimensions
%               are [3 3 nSym].
% symInfo       Structure containing additional information about the space
%               group with the following fields:
%   name            Name of the space group in string. If function called
%                   with no input, name stores the name of all spase groups
%                   from symmetry.dat in a cell.
%   str             The string of the symmetry operations.
%   num             The index of the symmetry in the symmetry.dat file.
%
%
% See also SW, SWSYM.POSITION, SPINW.ATOM, SPINW.MATOM, SWSYM.GENERATOR,
% SWSYM.POINT.
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