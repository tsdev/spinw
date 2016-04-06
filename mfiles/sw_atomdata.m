function varargout = sw_atomdata(atomSymb, datType)
% returns information on elements stored in the atom.dat file
%
% [data atomLabel] = SW_ATOMDATA(atomSymb, dataType)
%
% Input:
%
% atomSymb  String of the name of the atom, for example 'He'. If the string
%           contains whitespace character, the second word will be used to
%           identify the atom.
% dataType  Type of information requested:
%               radius      Atomic radius.
%               color       Color of the atom from the CPK color scheme.
%               mass        Average mass of the element.
%               longname    Name of the element.
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

if nargin == 0
    help sw_atomdata
    return
end

% read in the atom definition file
rPath = [sw_rootdir 'dat_files' filesep 'atom.dat'];
atom  = sw_readtable(rPath);

% split multiple words and use the second word if exists
atomSymb = strword(atomSymb,2,true);
atomSymb = atomSymb{1};

% cut M from the beginning of the atom label
upStr = isstrprop(atomSymb,'upper');
if (numel(atomSymb)>=2) && all(upStr(1:2))
    atomSymb = atomSymb(2:end);
end
% cut end symbols
atomSymb = atomSymb(isstrprop(atomSymb,'alpha'));

% find atom symbol
idx = find(strcmp({atom(:).name},atomSymb));

if isempty(idx)
    atom       = struct;
    atom.name  = '';
    atom.R     = 1;
    atom.RGB   = [255 167 0];
    atom.mass  = 0;
    atom.longname = '';
    
else
    atom = atom(idx);
end

switch datType
    case 'radius'
        data = atom.R;
    case 'color'
        data = atom.RGB;
    case 'mass'
        data = atom.mass;
    case 'longname'
        data = atom.longname;
    otherwise
        error('sw_atomdata:WrongInput','datType has to be one of the string options, see help sw_atomdata!');
end

varargout{1} = data;

if nargout > 1
    varargout{2} = atomSymb;
end

end