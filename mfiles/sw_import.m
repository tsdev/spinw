function obj = sw_import(fName, toPlot)
% create SpinW object from FullProf Studio file (.fst)
%
% At present can only import and plot crystal structure data.
%
% obj = SW_IMPORT(fName, {toPlot})
%
% Input:
%
% fName     String, contain the file location and name of the .fst file.
% toPlot    If true the structure will be plotted, default is false.
%

if nargin == 0
    help sw_import
    return
end

if nargin < 2
    toPlot = false;
end

% read the file into a string
str = fileread(fName);
% split string up into lines

str = strsplit(str,'\n');
% remove empty lines
str = str(~cellfun(@(C)isempty(C),str));
% store all imported data
dat = struct('atom',struct('label',{},'r',{}),'cell',[],'box',[],'matom',struct('label',{},'r',{},'M',{}),'k',zeros(3,0));

for ii = 1:numel(str)
    sstr = str{ii};
    
    if ismember(sstr(1),'!{}')
        % !COMMENT and grouping
        continue
    end
    
    sSplit = strsplit(sstr,' ');
    switch sSplit{1}
        case 'SPACEG'
            dat.spgr = sstr(8:end);
        case 'CELL'
            dat.cell = str2double(sSplit(2:7));
        case 'BOX'
            dat.box = str2double(sSplit(2:7));
        case 'ATOM'
            dat.atom(end+1).label = [sSplit{2} ' ' sSplit{3}];
            dat.atom(end).r       = str2double(sSplit(4:6))';
        case 'K'
            dat.k = [dat.k str2double(sSplit(2:4))'];
        case 'BKG'
            % [RGBT]
            dat.bkg = str2double(sSplit(2:5))*255;
        case 'MATOM'
            if numel(sSplit{3})>1
                sSplit{3}(2) = lower(sSplit{3}(2));
            end
            dat.matom(end+1).label = [sSplit{2} ' ' sSplit{3}];
            dat.matom(end).r       = str2double(sSplit(4:6))';
        case 'SKP'
            % M: 3 x 1 x nK
            kIdx = str2double(sSplit{2});
            Mk   = 0.5*(str2double(sSplit(4:6))+1i*str2double(sSplit(7:9)))*exp(-2*pi*1i*str2double(sSplit{10}));
            dat.matom(end).M(:,1,kIdx) = transpose(Mk);
    end
end

% create the SpinW model
obj0 = spinw;
obj0.genlattice('lat_const',dat.cell(1:3),'angled',dat.cell(4:6),'spgr',dat.spgr);
if ~isempty(dat.atom)
    obj0.addatom('r',[dat.atom(:).r],'label',{dat.atom(:).label});
elseif ~isempty(dat.matom)
    % add fictious spin to the atom
    obj0.addatom('r',[dat.matom(:).r],'label',{dat.matom(:).label},'S',ones(1,numel(dat.matom)));
end

if ~isempty(dat.k)
    % generate the single-k magnetic structure
    obj0.mag_str.k = dat.k(:,1)';
    % determine normal vector
    obj0.mag_str.n = cross(real(dat.matom(1).M),imag((dat.matom(1).M)))';
    obj0.mag_str.n = obj0.mag_str.n./norm(obj0.mag_str.n);
    
    % determine real magnetic moment components
    M0 = [dat.matom(:).M];
    % number of k-vectors
    nk = size(dat.k,2);
    % convert from the lattice component coordinate system to Descartes
    % coordinate system
    for ii = 1:nk
        M0(:,:,ii) = obj0.basisvector(1)*M0(:,:,ii);
    end
    
    % keep only the first wave vector
    % the given magnetic moments are in Bohr magneton, SpinW assumes g=2, thus
    % the moments are divided by this number
    g0 = 2;
    % calculate the real part of the magnetic structure
    obj0.mag_str.S = 2*real(M0(:,:,1))/g0;
end

if nargout > 0
    obj = obj0;
end

if toPlot
    plot(obj0,'range',reshape(dat.box',[2 3])')
end

end