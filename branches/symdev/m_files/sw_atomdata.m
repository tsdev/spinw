function data = sw_atomdata(atomSymb, datType)
% data = SW_ATOMDATA(atomSymb, dataType) returns information about atoms
% from the atom.dat file.
%
% Input:
%
% atomSymb  String of the name of the atom, for example 'He'.
% datType   Type of information requested:
%               radius  Atomic radius.
%               color   Color of the atom from the CPK color scheme.
%
% Example:
% sw_atomdata('H','radius') = 0.37
% If the atom label does not exists, the function returns radius = 1, 
% color = [255 167 0].
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

% cut M from the beginning of the atom label
upStr = isstrprop(atomSymb,'upper');
if (numel(atomSymb)>=2) && all(upStr(1:2))
    atomSymb = atomSymb(2:end);
end
% cut end symbols
atomSymb = atomSymb(isstrprop(atomSymb,'alpha'));

atom = struct;
idx = 1;
found = false;
while (~feof(fid)) && ~found
    fLine = fgets(fid);
    [atom.name, ~, ~, nextIdx] = sscanf(fLine,'%s',1);
    [atom.rad, ~, ~, nextIdx2] = sscanf(fLine(nextIdx:end),'%f,',1);
    atom.col = sscanf(fLine((nextIdx+nextIdx2):end),'%f ',[1 3]);
    found = strcmpi(atom.name,atomSymb);
    idx = idx + 1;
end

fclose(fid);

if strcmpi('radius',datType)
    if found
        data = atom.rad;
    else
        data = 1;
    end
elseif strcmpi('color',datType)
    if ~isempty(atom.col)
        data = atom.col;
    else
        data = [255 167 0];
    end
else
    error('sw_atomdata:WrongInput','datType has to be either radius or color!');
end

end