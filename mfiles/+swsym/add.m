function sym = add(symStr, symName)
% saves user defined symmetry operators
%
% sym = SWSYM.ADD(symStr, {symName})
%
% It saves the symmetry generators in symStr into the symmetry.dat file and
% returns the line number of the space group in the symmetry.dat file.
%
% Input:
% symStr        String, that contains the operators of the space group. If
%               not only the generators are given, a possible set of
%               generators will be determined and only those will be saved.
% symName       Label for the space group.
%
% See also SW_GENSYM, SW_SYMGETGEN.
%

if nargin == 0
    help sw_addsym
    return
end

symPath = [sw_rootdir 'dat_files' filesep 'symmetry.dat'];
% Count the number of lines
fid = fopen(symPath);
if fid == -1
    error('add:FileNotFound',['Symmetry definition file not found: '...
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
symOp    = swsym.generator(symStr);
[~, isG] = swsym.genreduce(symOp);

% parse the input string
idx = 1;
parStr = {};
while ~isempty(symStr)
    [parStr{idx}, symStr] = strtok(symStr,';'); %#ok<STTOK,AGROW>
    idx = idx + 1;
end
parStr  = parStr(isG);
parStr  = parStr(:)';
parStr  = [parStr; repmat({';'},1,numel(parStr))];
symStrG = [parStr{:}];
symStrG = symStrG(1:end-1);

% Write the file
fid = fopen(symPath,'a');
fprintf(fid,'\n%4d  %-11s: %s',nLines+1,symName,symStrG);
fclose(fid);

sym = nLines + 1;

end