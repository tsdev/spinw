function varargout = sw_atomdata(atomSymb, datType)
% returns information on elements stored in the atom.dat file
%
% [data atomLabel] = SW_ATOMDATA(atomSymb, {dataType})
%
% Input:
%
% atomSymb  String of the name of the atom, for example 'He' or the atomic
%           number Z. If the string contains whitespace character, the
%           second word will be used to identify the atom.
% dataType  Type of information requested:
%               name        Atomic symbol.
%               radius      Atomic radius.
%               color       Color of the atom from the CPK color scheme.
%               mass        Average mass of the element.
%               longname    Name of the element.
%               Z           tomic index.
%               all         All atomic data returned in a struct.
%
% Example:
% sw_atomdata('H','radius') = 0.37
% If the atom label does not exists, the function returns radius = 1,
% color = [255 167 0].
%
% optional second output is 'atomLabel' that contains the name of the atom
% clean.
%
% See also SW_MFF.
%

if nargin == 1
    datType = 'all';
end

if nargin == 0
    help sw_atomdata
    return
end

% read in the atom definition file
rPath = [sw_rootdir 'dat_files' filesep 'atom.dat'];
atom  = sw_readtable(rPath);

% add atomic number for the first 111 atom (original definitions)
cellIdx = num2cell(1:111);
% FANCY
[atom(1:111).Z] = deal(cellIdx{:});

% add empty placeholder atom
atom(end+1).name   = '';
atom(end).R        = 1;
atom(end).RGB      = [255 167 0];
atom(end).mass     = 0;
atom(end).longname = '';
atom(end).Z        = 0;

if ischar(atomSymb)
    atomSymb = {atomSymb};
end

if iscell(atomSymb)
    % point to the empty atom data by default
    idx = zeros(1,numel(atomSymb))+numel(atom);
    
    for ii = 1:numel(atomSymb)
        % split multiple words and use the second word if exists
        atomSymb0 = strword(atomSymb{ii},2,true);
        atomSymb0 = atomSymb0{1};
        
        % cut M from the beginning of the atom label
        upStr = isstrprop(atomSymb0,'upper');
        if (numel(atomSymb0)>=2) && all(upStr(1:2))
            atomSymb0 = atomSymb0(2:end);
        end
        % cut end symbols
        atomSymb0 = atomSymb0(isstrprop(atomSymb0,'alpha'));
        % find atom symbol
        idx0 = find(strcmp({atom(:).name},atomSymb0));
        
        atomSymb{ii} = atomSymb0;
        
        if ~isempty(idx0)
            idx(ii) = idx0;
        end
    end
    
    if numel(atomSymb) == 1
        atomSymb = atomSymb{1};
    end
else
    idx = atomSymb;
end

% select atoms
atom = atom(idx);

switch datType
    case 'name'
        data = {atom.name};
    case 'radius'
        data = [atom.R];
    case 'color'
        data = reshape([atom.RGB],3,[]);
    case 'mass'
        data = [atom.mass];
    case 'longname'
        data = {atom.longname};
    case 'Z'
        data = [atom.Z];
    case 'all'
        data = rmfield(atom,'MODE');
    otherwise
        error('sw_atomdata:WrongInput','datType has to be one of the string options, see help sw_atomdata!');
end

varargout{1} = data;

if nargout > 1
    varargout{2} = atomSymb;
end

end