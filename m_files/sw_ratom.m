function radius = sw_ratom(atomSymb)
% radius = SW_RATOM(atomSymb) returns the atomic radius based on atom
% symbol. Atomic radius data is stored in atomradius.dat file.
%
% Example:
% sw_ratom('H') = 0.37
% If the atom label does not exists, the function returns -1.
%
% See also SW_MFF.
%

% Open the atomic radius definition file.
fid = fopen([sw_rootdir 'atomradius.dat']);
if fid == -1
    
    error('spinw:sw_ratom:FileNotFound',['Atomic radius definition file not found: '...
        regexprep(sw_rootdir,'\' , '\\\') 'atomradius.dat!']);
end

atom = struct;
idx = 1;
found = false;
while (~feof(fid)) && ~found
    fLine = fgets(fid);
    [atom.name, ~, ~, nextIdx] = sscanf(fLine,'%s',1);
    atom.rad = sscanf(fLine(nextIdx:end),'%f,',1);
    found = strcmpi(atom.name,atomSymb);
    idx = idx + 1;
end

fclose(fid);

if found
    radius = atom.rad;
else
    radius = -1;
end

end