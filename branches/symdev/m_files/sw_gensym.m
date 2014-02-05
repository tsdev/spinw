function [symOp, symTr, symName, symStr, symNum] = sw_gensym(varargin)
% [symOp, symTr, symName, symStr, symNum] = SW_GENSYM({sym}) gives the
% symmetry elements based on the space group number or given list of
% symmetry operators. Without arguments, returns the name of all space
% groups stored in symmetry.dat.
%
% Input:
%
% sym           It is either the name of the space group or the index from
%               the International Tables of Crystallography.
%
% Output:
%
% symOp         Rotation matrices, dimensions are [3 3 nOp].
% symTr         Translation vectors, dimensions are [3 nOp].
% symName       Name of the space group, stored in cells.
% symStr        The string of the symmetry operations.
% symNum        The index of the symmetry in the symmetry.dat file.
%
% See also SW_ADDSYM, SW, SW.GENCOUPLING, SW_GENCOORD.
%

if nargin == 0
    help sw_gensym;
    return;
end

% Open the symmetry definition file.
symPath = [sw_rootdir 'dat_files' filesep 'symmetry.dat'];
fid = fopen(symPath);
if fid == -1
    
    error('spinw:sw_gensym:FileNotFound',['Symmetry definition file not found: '...
        regexprep(symPath,'\' , '\\\') '!']);
end

% Just returns the name of all space groups.
if nargin == 0
    ii = 1;
    symName = {};
    while ~feof(fid)
        textLine    = fgetl(fid);
        symName{ii} = [textLine(7:17) sprintf(' (%3i)',ii)]; %#ok<AGROW>
        ii = ii+1;
    end
    symOp = [];
    symTr = [];
    symStr  = [];
    fclose(fid);
    return

elseif nargin == 1
    
    if iscell(varargin{1})
        varargin{1} = varargin{1}{1};
    end
    
    if ischar(varargin{1})
        if isempty(varargin{1})
            symStr = 'x,y,z';
            symName = 'P 1';
        else
            % find symmetry label
            symName = varargin{1};
            symName(end+1:11) = 32;
            symIdx = 0;
            ii     = 1;
            while (symIdx == 0) && ~feof(fid)
                textLine    = fgetl(fid);
                if strfind(symName,textLine(7:17))
                    symIdx = ii;
                end
                ii = ii+1;
            end
            if symIdx == 0
                fclose(fid);
                error('sw:sw_gensym:WrongInput','Symmetry name does not exists (case insensitive)!');
            end
            symNum = symIdx;
            symStr     = textLine(20:end);
        end
        
    else
        symNum = varargin{1};
        if symNum<=0
            fclose(fid);
            error('spinw:sw_gensym:WrongInput','Symmetry number has to be positive integer!');
        end
        ii = 1;
        while (ii<=symNum) && ~feof(fid)
            textLine = fgetl(fid);
            ii = ii+1;
        end
        fclose(fid);
        if ii <= symNum
            fclose(fid);
            error('spinw:sw_gensym:WrongInput','Symmetry number not found!')
        end
        symStr   = textLine(20:end);
        symName = textLine(7:17);
    end
elseif nargin >= 2
    symName = varargin{1};
    symStr   = varargin{2};

end

symOp = zeros(3,3,30);
symTr = zeros(3,30);
vNew   = zeros(3,1);

nNew  = 1;
nOp   = 1;
nSign = 1;

ii=1;
while(ii<=length(symStr))
    if symStr(ii)==','
        symOp(nNew,:,nOp) = vNew;
        vNew  = vNew*0;
        nSign = 1;
        nNew  = mod(nNew,3)+1;
    elseif symStr(ii)==';'
        symOp(nNew,:,nOp) = vNew;
        vNew = vNew*0;
        nSign = 1;
        nNew  = 1;
        nOp   = nOp+1;
    elseif symStr(ii)=='x'
        vNew(1) = nSign;
    elseif symStr(ii)=='y'
        vNew(2) = nSign;
    elseif symStr(ii)=='z'
        vNew(3) = nSign;
    elseif symStr(ii)=='-'
        nSign = -1;
    elseif symStr(ii)=='+'
        nSign = 1;
    elseif (symStr(ii)=='1')||(symStr(ii)=='2')||(symStr(ii)=='3')
        symTr(nNew,nOp) = (symStr(ii)-'0')/(symStr(ii+2)-'0');
        ii = ii+2;
    end
    ii = ii+1;
    
end

symOp(nNew,:,nOp) = vNew;
symOp = symOp(:,:,1:nOp);
symTr = symTr(:,1:nOp);

% cut trailing spaces from symName
if isnan(symName)
    symName = '';
else
    symName = symName(1:find(diff([symName '  ']==32),1,'last'));
end

    

end