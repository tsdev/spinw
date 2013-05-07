function [transf, transl, symName, symOp, symNumber] = sw_gensym(varargin)
% [transf, transl, symName, symOp] = SW_GENSYM({sym}) gives the symmetry elements
% based on the space group number or given list of symmetry operators.
% Without arguments, returns the name of all space groups stored in
% symmetry.dat.
%
% Input:
%
% sym           It is either the name of the space group or the index from
%               the International Tables of Crystallography.
%
% Output:
%
% tranf         Rotation matrices, dimensions are [3 3 nOp].
% tranl         Translation vectors, dimensions are [3 nOp].
% symName       Name of the space group, stored in cells.
% symOp         The string of the symmetry operations.
% symNumber     The index of the symmetry in the symmetry.dat file.
%
% See also SW_ADDSYM, SW, SW.GENCOUPLING.
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
    transf = [];
    transl = [];
    symOp  = [];
    fclose(fid);
    return
end

if iscell(varargin{1})
    varargin{1} = varargin{1}{1};
end

if ischar(varargin{1})
    if isempty(varargin{1})
        symOp = 'x,y,z';
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
            error('sw:sw_gensym:WrongInput','Symmetry name does not exists (case insensitive)!');
        end
        symNumber = symIdx;
        symOp     = textLine(20:end);
    end
    
else
    symNumber = varargin{1};
    if symNumber<=0
        fclose(fid);
        error('spinw:sw_gensym:WrongInput','Symmetry number has to be positive integer!');
    end
    ii = 1;
    while (ii<=symNumber) && ~feof(fid)
        textLine = fgetl(fid);
        ii = ii+1;
    end
    fclose(fid);
    if ii <= symNumber
        error('spinw:sw_gensym:WrongInput','Symmetry number not found!')
    end
    symOp   = textLine(20:end);
    symName = textLine(7:17);
end

transf = zeros(3,3,10);
transl = zeros(3,10);
vNew   = zeros(3,1);

nNew  = 1;
nOp   = 1;
nSign = 1;

ii=1;
while(ii<=length(symOp))
    if symOp(ii)==','
        transf(nNew,:,nOp) = vNew;
        vNew  = vNew*0;
        nSign = 1;
        nNew  = mod(nNew,3)+1;
    elseif symOp(ii)==';'
        transf(nNew,:,nOp) = vNew;
        vNew = vNew*0;
        nSign = 1;
        nNew  = 1;
        nOp   = nOp+1;
    elseif symOp(ii)=='x'
        vNew(1) = nSign;
    elseif symOp(ii)=='y'
        vNew(2) = nSign;
    elseif symOp(ii)=='z'
        vNew(3) = nSign;
    elseif symOp(ii)=='-'
        nSign = -1;
    elseif symOp(ii)=='+'
        nSign = 1;
    elseif (symOp(ii)=='1')||(symOp(ii)=='2')||(symOp(ii)=='3')
        transl(nNew,nOp) = (symOp(ii)-'0')/(symOp(ii+2)-'0');
        ii = ii+2;
    end
    ii = ii+1;
    
end

transf(nNew,:,nOp) = vNew;
transf = transf(:,:,1:nOp);
transl = transl(:,1:nOp);

% cut trailing spaces from symName
symName = symName(1:find(diff([symName '  ']==32),1,'last'));

end