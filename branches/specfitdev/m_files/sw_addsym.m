function sym = sw_addsym(symStr, symName)
% sym = SW_ADDSYM(symStr, {symName}) saves the symmetry generators in
% symStr into the symmetry.dat file and returns the number of the space
% group.
%
% See also SW_GENSYM.
%

symPath = [sw_rootdir 'dat_files' filesep 'symmetry.dat'];
% Count the number of lines
fid = fopen(symPath);
if fid == -1
    error('spinw:sw_gensym:FileNotFound',['Symmetry definition file not found: '...
        regexprep(symPath,'\' , '\\\') '!']);
end

nLines = 0;
while (fgets(fid) ~= -1)
    nLines = nLines+1;
end
fclose(fid);

if nargin == 1
    symName = ['sym' num2str(nLines+1)];
end

% Write the file
fid = fopen(symPath,'a');
fprintf(fid,'\n%4d  %-11s: %s',nLines+1,symName,symStr);
fclose(fid);

sym = nLines + 1;
end