function varargout = sw_atomdata(atomSymb, datType)
% [data atomLabel] = SW_ATOMDATA(atomSymb, dataType) returns information about atoms
% from the atom.dat file.
%
% Input:
%
% atomSymb  String of the name of the atom, for example 'He'. If the string
%           contains whitespace character, the second word will be used to
%           identify the atom.
% datType   Type of information requested:
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
    help sw_atomdata;
    return;
end

% Open the atomic radius definition file.
rPath = [sw_rootdir 'dat_files' filesep 'atom.dat'];
fid = fopen(rPath);
if fid == -1
    error('spinw:sw_atomdata:FileNotFound',['Atom definition file cannot be found: '...
        regexprep(rPath,'\' , '\\\') '!']);
end

% split multiple words and use the second word if exists
atomSymb = strword(atomSymb,2,true);
atomSymb = atomSymb{1};

% atomSymb = strsplit(atomSymb);
% wordIdx = find(cellfun(@numel,atomSymb));
% if numel(wordIdx) > 1
%     atomSymb = atomSymb{wordIdx(2)};
% else
%     atomSymb = atomSymb{wordIdx(1)};
% end

% cut M from the beginning of the atom label
upStr = isstrprop(atomSymb,'upper');
if (numel(atomSymb)>=2) && all(upStr(1:2))
    atomSymb = atomSymb(2:end);
end
% cut end symbols
atomSymb = atomSymb(isstrprop(atomSymb,'alpha'));

atom = struct;

% read in all data
% Symbol IonicRadius R G B Mass Name
aData = textscan(fid,'%s %f %d %d %d %f %s');
fclose(fid);

% find atom symbol
idx = find(strcmp(aData{1},atomSymb));

if isempty(idx)
    found = false;
    atom.name  = '';
    atom.rad   = 1;
    atom.color = [255 167 0];
    atom.mass  = 0;
    atom.longname = '';
    
else
    found = true;
    atom.name  = aData{1}{idx};
    atom.rad   = aData{2}(idx);
    atom.color = double([aData{3}(idx) aData{4}(idx) aData{5}(idx)]);
    atom.mass  = aData{6}(idx);
    atom.longname = aData{7}{idx};
end

switch datType
    case 'radius'
        data = atom.rad;
    case 'color'
        data = atom.color;
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