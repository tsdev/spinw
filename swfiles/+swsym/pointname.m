function [intStr, info] = pointname(O,r)
% determine the name of the point group
%
% symStr = SWSYM.POINTNAME(symOp,{r})
%
% Input:
%
% symOp     Stack of symmetry operators in a matrix with dimensions of
%           [3,4,nOp].
% r         Position in the unit cell, if symOp contains translation
%           vectors. Default value is [0;0;0].
%
% Output:
%
% symStr    Name of the point group.
%

if nargin<2
    r = [0;0;0];
end

% generate the point group operators corresponding to the given position
symOp = swsym.point(O,r);

% remove unit operator
[~,symSel] = unique(reshape(cat(3,eye(3),symOp),9,[])','rows');
symOp = symOp(:,:,symSel(symSel>1)-1);

% remove higher order operators
idx = 1;
while idx < size(symOp,3)
    oGen = swsym.oporder(symOp(:,:,idx));
    if oGen > 2
        sPow = 2:(oGen-1);
        % generate all powers
        nGen = oGen-2;
        symGen = zeros(3,3,nGen);
        for ii = sPow
            symGen(:,:,ii-1) = symOp(:,:,idx)^ii;
        end
        % calculate unique matrices
        [~,symSel] = unique(reshape(cat(3,symGen,symOp),9,[])','rows');
        symOp = symOp(:,:,symSel(symSel>nGen)-nGen);
    end
    idx = idx + 1;
end

% determine the symmetry operator names


% get the symmetry operator names
if isempty(symOp)
    info = struct('name','1','name2','E','axis',zeros(3,1),'order',1,'isrotation',false);
else
    info = swsym.pointopname(symOp);
end

% load the point group name database
pDat = sw_readtable([sw_rootdir 'dat_files' filesep 'pointgroup.dat']);

% The Hermann?Mauguin or international notation:
% n                 The principal axis (z axis) is a rotation axis of order n
% -n                The principal axis (z axis) is a rotoinversion axis of order n.
% n2 or -n2         A binary axis perpendicular to the principal axis.
% nm or -nm         A mirror (reflection plane) parallel to the principal axis.
% n/m or -n/m       A mirror perpendicular to the principal axis.
% Further entries   They refer to secondary axes.


% string storing the space group name
intStr = '';
idx = 1;
% number of point group operators
nPoint = numel(info);
% check first the principal axis, then 2 other secondary axes
while idx<=nPoint && info(idx).isrotation && idx<4
    intStr = [intStr info(idx).name];
    if info(idx).isrotation && idx<nPoint && info(idx).order>2
        % is there a binary axis perpendicular to the principal axis?
        % use 1/2 for empty entries
        axProd = [ones(1,idx)/2 abs(sum(info(idx).axis.*[info((idx+1):end).axis],1))];
        if any(ismember({info(axProd<10*eps).name},'2'))
            intStr = [intStr '2'];

        % is there any mirror plane parallel to the principal axis?
        elseif any(ismember({info(axProd<10*eps).name},'m'))
            intStr = [intStr 'm'];
        
        % is there any mirror plane perpendicular to the principal axis?
        elseif any(ismember({info(axProd>(1-10*eps)).name},'m'))
            intStr = [intStr '/m'];
        end
    end
    idx = idx+1;
    intStr = [intStr ' '];
end

if isempty(intStr)
    intStr = info(1).name;
end

return

multi = zeros(1,10);
multi(info.type) = info.multi;

% type name    trace     order
%  1 E     --> chi =  3, N = 1
%  2 i     --> chi = -3, N = 2
%  3 sigma --> chi =  1, N = 2   == S2
%  4 C2    --> chi = -1, N = 2
%  5 C3    --> chi =  0, N = 3
%  6 C4    --> chi =  1, N = 4
%  7 C6    --> chi =  2, N = 6
%  8 S3    --> chi = -2, N = 6
%  9 S4    --> chi = -1, N = 4
% 10 iC3   --> chi =  0, N = 6

% missing:
% S4	-4
% C3i	-3 --> S6
% T	23
% Th	m-3
% O	432

pStr = [];

% determines the point group

if any(info.multi(ismember(info.type,5:7))>1)
    % Does if have 2 or more Cn (n>2)?
    if any(info.type==2)
        % Does it have center of inversion?
        pStr = {'Td' '-43m'};
    else
        pStr = {'Oh' 'm-3m'};
    end
else
    % Does it have Cn axis?
    if any(ismember(info.type,4:7))
        % Are there n C2 axis perpendicular to the principal?
        n = max(info.type(ismember(info.type,4:7)));
        ns = num2str(n);
        if (n==4 && info.multi(info.type==4)==3) || (n==7 && info.multi(info.type==4)==6) ||...
                (info.multi(info.type==4)==(n-2))
            
            if multi(3) == 1
                % Is there a horizontal mirror plane?
                pStr = ['D' ns 'h'];
                switch n
                    case 2
                        pStr{2} = 'mmm';
                    case 3
                        pStr{2} = '-62m';
                    case 4
                        pStr{2} = '4/mmm';
                    case 6
                        pStr{2} = '6/mmm';
                end
            elseif multi(3) == n
                % Are there n dihedral mirror plane?
                pStr = ['D' ns 'd'];
                switch n
                    case 2
                        pStr{2} = '-42m';
                    case 3
                        pStr{2} = '-3m';
                end
                
            else
                pStr = ['D' ns];
                switch n
                    case 2
                        pStr{2} = '222';
                    case 3
                        pStr{2} = '32';
                    case 4
                        pStr{2} = '422';
                    case 6
                        pStr{2} = '622';
                end
            end
        elseif multi(3) > 0
            % Is there a horizontal mirror plane?
            pStr = 'Cnh';
            switch n
                case 2
                    pStr{2} = '2/m';
                case 3
                    pStr{2} = '-6';
                case 4
                    pStr{2} = '4/m';
                case 6
                    pStr{2} = '6/m';
            end
        elseif multi(3)==n
            % Are there n vertical mirror planes? -> Cnv
            pStr = 'Cnv';
            switch n
                case 2
                    pStr{2} = 'mm2';
                case 3
                    pStr{2} = '3m';
                case 4
                    pStr{2} = '4mm';
                case 6
                    pStr{2} = '6mm';
            end
        elseif multi(9)>0
            % Is  there an S2n axis?
            pStr = 'S2n';
        else
            pStr = {['C' ns] ns};
        end
    else
        if any(ismember(info.type,3))
            % Does it have a mirror plane?
            pStr = {'Cs' 'm'};
        elseif any(ismember(info.type,2))
            % Does it have a center of inversion?
            pStr = {'Ci' '-1'};
        else
            pStr = {'C1' '1'};
        end
    end
end


end