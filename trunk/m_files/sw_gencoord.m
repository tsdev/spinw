function [symOp, symTr] = sw_gencoord(sym, fid, tol)
% [symOp, symTr] = SW_GENCOORD(sym, {fid}, {tol}) calculates the general
% coordinates for a given space group.
%
% Input:
%
% sym           Line index in the symmetry.dat file or string of the
%               symmetry operators.
% {fid}         For printing the symmetry operators:
% {print}       Optional input, if true, the result will be plotted on the
%               command window.
%
% Output:
%
% symOp         The rotational part of the symmetry operators, dimensions
%               are [3 3 nSym]. 
% symTr       	The translation part of the symmetry operators, dimensions
%               are [3 nSym].
%
% See also SW, SW.ATOM, SW.MATOM, SW_GENCOUPLING.
%

if nargin == 0
    help sw_gencoord;
    return;
end

if nargin < 3
    tol = 1e-5;
end
if nargin < 2
    fid = 0;
end

[genOp, transl, symName] = sw_gensym(sym);

nGen  = size(genOp,3);
symOp = eye(3);
symTr = zeros(3,1);

% determine the order of the generators
N = zeros(1,nGen);
M = zeros(1,nGen);
P = zeros(1,nGen);

% generate all symmetry elements
for ii = 1:nGen
    % order of the rotational matrix
    N(ii) = sw_rotorder(genOp(:,:,ii));
    % order of the translation
    M(ii) = max(12./gcd(abs(transl(:,ii))*12,12));
    if M(ii) == 0
        M(ii) = 1;
    end
    % take the least common multiplier
    P(ii) = lcm(M(ii),N(ii));
    for jj = 1:(P(ii)-1)
        R = genOp(:,:,ii)^jj;clc
        
        nSym = size(symOp,3);
        for kk = 1:nSym
            RS = symOp(:,:,kk)*R;
            idx = permute(sum(sum(abs(bsxfun(@minus,symOp,RS)),1),2),[3 1 2]) > tol;
            % adds new operator to the list if it differs from all
            if all(idx)
                symOp = cat(3,symOp,RS);
                symTr = [symTr mod(transl(:,ii)*jj+symTr(:,kk),1)]; %#ok<AGROW>
            end
            
        end
    end
end

% print symmetry elements
if fid~=0
    fprintf(fid, 'General coordinates of space group: %s\n',symName);
    
    nSym = size(symOp,3);
    for ii = 1:nSym
        fprintf(fid,'(%d) ',ii);
        for jj = 1:3
            first = true;
            for kk = 1:3
                symE = symOp(jj,kk,ii);
                switch sign(symE)
                    case 1
                        if first
                            fprintf(fid,char(119+kk));
                        else
                            fprintf(fid,['+' 119+kk]);
                            first = false;
                        end
                    case -1
                        fprintf(fid,['-' 119+kk]);
                        first = false;
                end
            end
            trE  = symTr(jj,ii);
            ndTr = [trE 1]*12/gcd(trE*12,12);
            if trE
                fprintf(fid,'%+d/%d',ndTr);
            end
            if jj < 3
                fprintf(fid,',');
            end
        end
        fprintf(fid,'\n');
    end
end

end