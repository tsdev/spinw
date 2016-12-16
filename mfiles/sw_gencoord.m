function [symOp, symTr, symName] = sw_gencoord(sym, fid, tol)
% calculates all symmetry operators for a given space group
%
% [symOp, symTr, symName] = SW_GENCOORD(sym, fid) 
%
% Input:
%
% sym           Line index in the symmetry.dat file or string of the
%               symmetry operators or cell: {symOp, symTr}.
%               For example:
%                   sym = 'P b n m';
%               Can be also a matrix with dimensions [3 4 nOp] in a
%               symmetry operator format.
% fid           For printing the symmetry operators:
%                   0   no printed output (Default)
%                   1   standard output (Command Line)
%                   fid text file opened before with the fid = fopen(path)
%
% Output:
%
% symOp         The rotational part of the symmetry operators, dimensions
%               are [3 3 nSym].
% symTr       	The translation part of the symmetry operators, dimensions
%               are [3 nSym].
% symName       String, the name of  the space group.
%
% See also SW, SW.ATOM, SW.MATOM, SW.GENCOUPLING, SW_POINTSYM, SW_GENATPOS,
% SW_GENSYM, SW_ISSYMOP.
%

if nargin == 0
    help sw_gencoord
    return
end

% tolerance for numerical error
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
P = zeros(1,nGen);

% generate all symmetry elements
for ii = 1:nGen
    R0 = genOp(:,:,ii);
    T0 = transl(:,ii);
    % order of the symmetry operator
    P(ii) = sw_symorder(R0,T0);
    
    R = eye(3);
    T = zeros(3,1);
    
    for jj = 1:(P(ii)-1)
        R = R0*R;
        T = R0*T + T0;
        
        nSym = size(symOp,3);
        for kk = 1:nSym
            RS = R*symOp(:,:,kk);
            %TS = mod(R*symTr(:,kk) + T,1);
            TS = mod(round(24*(R*symTr(:,kk)+T))/24,1);
            
            idxR = permute(sum(sum(abs(bsxfun(@minus,symOp,RS)),1),2),[3 1 2]) > tol;
            idxT = sum(abs(bsxfun(@minus,symTr,TS)),1)' > tol;
            
            % adds new operator to the list if it differs from all
            if all(idxR | idxT)
                symOp = cat(3,symOp,RS);
                symTr = [symTr TS]; %#ok<AGROW>
            end
            
        end
    end
end

% print symmetry elements
if fid~=0
    fprintf(fid, 'General coordinates of space group: %s\n',symName);
    
    sStr = strtrim(strsplit(sw_gensymstr(symOp,symTr),';')');
    for ii = 1:numel(sStr)
        fprintf(fid,'(%02d) %s\n',ii,sStr{ii});
    end
    
%     nSym = size(symOp,3);
%     for ii = 1:nSym
%         fprintf(fid,'(%d) ',ii);
%         for jj = 1:3
%             first = true;
%             for kk = 1:3
%                 symE = symOp(jj,kk,ii);
%                 switch sign(symE)
%                     case 1
%                         if first
%                             fprintf(fid,char(119+kk));
%                         else
%                             fprintf(fid,['+' 119+kk]);
%                             first = false;
%                         end
%                     case -1
%                         fprintf(fid,['-' 119+kk]);
%                         first = false;
%                 end
%             end
%             trE  = symTr(jj,ii);
%             ndTr = [trE 1]*12/gcd(trE*12,12);
%             if trE
%                 fprintf(fid,'%+d/%d',ndTr);
%             end
%             if jj < 3
%                 fprintf(fid,',');
%             end
%         end
%         fprintf(fid,'\n');
%     end
end

end