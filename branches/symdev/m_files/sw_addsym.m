function sym = sw_addsym(symStr, symName)
% sym = SW_ADDSYM(symStr, {symName}) saves the symmetry generators in
% symStr into the symmetry.dat file and returns the number of the space
% group.
%
% Input:
% symStr        Name for the 
% See also SW_GENSYM.
%

if nargin == 0
    help sw_addsym;
    return;
end

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

% determine the symmetry generators
[symOp, symTr] = sw_gensym(symName, symStr);
[~, ~, isG] = sw_symgetgen(symOp, symTr);
% parse the input string
idx = 1;
parStr = {};
while ~isempty(symStr)
    [parStr{idx}, symStr] = strtok(symStr,';'); %#ok<STTOK,AGROW>
    idx = idx + 1;
end
parStr = parStr(isG);
parStr = parStr(:)';
parStr = [parStr; repmat({';'},1,numel(parStr))];
symStrG = [parStr{:}];
symStrG = symStrG(1:end-1);

% Write the file
fid = fopen(symPath,'a');
fprintf(fid,'\n%4d  %-11s: %s',nLines+1,symName,symStrG);
fclose(fid);

sym = nLines + 1;
end