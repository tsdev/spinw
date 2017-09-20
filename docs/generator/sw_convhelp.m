function str = sw_convhelp(fName, foredit)
% try to convert old help to new MD format

if nargin<2
    foredit = true;
end

% remove
if any(fName==newline)
    %str   = [{''} strsplit(fName,newline)];
    str   = strsplit(fName,newline);
else
    str   = strsplit(help(fName),newline);
end

% if there is no H1 line add an empty one
if isempty(strtrim(str{1}))
    str = [{' '} str];
end

% remove leading spaces
str   = sw_rmspace(str);
% exchange common symbols
sText = {'Angstrom'     'hbar'     'alpha'     'beta'     'gamma'     'degree'     '\^-1'   'default is'};
cText = {'\\\\Angstrom' '\\\\hbar' '\\\\alpha' '\\\\beta' '\\\\gamma' '\\\\degree' '$\^{-1}$' 'default value is'};

for ii = 1:numel(sText)
    str = regexprep(str,sText{ii},cText{ii});
end


% create trimmed lines
strTr = strtrim(str);
% find sections
oldSections = {'^Example[s]?[:]?$' '^Input[s]?[:]?$' '^Option[s]?[:]?$' '^Output[s]?[:]?$'  '^See also'};
newSections = {'Examples' 'Input Arguments' 'Name-Value Pair Arguments'  'Output Arguments' 'See Also'};

sIdx     = zeros(2,numel(oldSections));

for ii = 1:numel(oldSections)
    tIdx = find(cellfun(@(C)~isempty(C),regexp(strTr,oldSections{ii},'once')));
    if isempty(tIdx)
        sIdx(1,ii) = 0;
    else
        sIdx(1,ii) = tIdx;
    end
    
end
sIdx(2,:) = 1:numel(oldSections);
sIdx = sortrows(sIdx')';
sIdx = [sIdx(:,logical(sIdx(1,:))) [numel(strTr)+1;0]];

% separate See also line
saIdx   = sIdx(1,sIdx(2,:) == 5);
if ~isempty(saIdx)
    seealso = str{saIdx};
    seealso = lower(strtrim(seealso(9:end)));
    if seealso(end)=='.'
        seealso = seealso(1:end-1);
    end
    seealso(seealso==',') = [];
    seealso = strsplit(seealso,' ');
    if numel(seealso)==1
        seealso1 = ['[' seealso{1} ']'];
    else
        seealso1 = sprintf('[%s] \\| ',seealso{1:end-1});
        seealso1 = [seealso1 '[' seealso{end} ']'];
    end
    % save back into string
    str{saIdx} = seealso1;
end


% find lines
lIdx = cell(1,(size(sIdx,2)-1));
sec = struct('title',cell(1,size(sIdx,2)-1),'text',[],'idx',[]);

for ii = 1:(size(sIdx,2)-1)
    if sIdx(2,ii) == 5
        % see also line
        lIdx{ii} = sIdx(1,ii):(sIdx(1,ii+1)-1);
    else
        lIdx{ii} = (sIdx(1,ii)+1):(sIdx(1,ii+1)-1);
    end
    sec(ii).title = newSections{sIdx(2,ii)};
    sec(ii).text  = str(lIdx{ii});
    sec(ii).text  = sec(ii).text(cellfun(@(C)~isempty(strtrim(C)),sec(ii).text));
    
    sec(ii).idx   = sIdx(2,ii);
end

% remove empty lines
[~,newidx] = sort([sec.idx]);
sec = sec(newidx);

% fix syntax line
if numel(str)>2
    syntax = ['`' lower(str{3}) '`'];
    isOpt  = strfind(syntax,'''option');
    if ~isempty(isOpt)
        syntax = [strtrim(syntax(1:isOpt-1)) 'Name,Value)`'];
    end
    % add syntax and description
    sec = [struct('title',{'Syntax' 'Description'},'text',{{syntax} str(5:(sIdx(1,1)-1))},'idx',{-2 -1}) sec];
else
    sec = [struct('title',{'Description'},'text',{str(5:(sIdx(1,1)-1))},'idx',{-1}) sec];
end


% fix obj line
if any([sec.idx]==2)
    objIdx = find(cellfun(@(C)~isempty(C),regexp(sec([sec.idx]==2).text,'^obj\s')),1);
    if ~isempty(objIdx)
        sec([sec.idx]==2).text{objIdx} = 'obj [spinw] object.';
    end
end

% convert arguments and options into definitions
for ii = 2:3
    if any([sec.idx] == ii)
        text0 = sec([sec.idx] == ii).text;
        % lines of definition
        defIdx   = cellfun(@(C)~isempty(C) && C(1)~=' ' ,text0);
        %spaceIdx = cellfun(@(C)find(C==' ',1),text0(defIdx));
        [spaceIdx, spaceEnd] = regexp(text0(defIdx),'\s+','once');
        spaceEnd = round(mean(cell2mat(spaceEnd)));
        text0    = [text0;cell(1,numel(text0))];
        fdefIdx  = find(defIdx);
        for jj = 1:numel(fdefIdx)
            text0{2,fdefIdx(jj)} = [':' text0{1,fdefIdx(jj)}(spaceEnd:end)];
            if foredit
                text0{1,fdefIdx(jj)} = [newline '% `' text0{1,fdefIdx(jj)}(1:(spaceIdx{jj}-1)) '`'];
            else
                text0{1,fdefIdx(jj)} = [newline '`' text0{1,fdefIdx(jj)}(1:(spaceIdx{jj}-1)) '`'];
            end
        end
        text0(1,~defIdx) = cellfun(@(C)[' ' C(spaceEnd:end)],text0(1,~defIdx),'UniformOutput',false);
        sec([sec.idx]==ii).text = text0(cellfun(@(C)~isempty(C),text0(:)));
        
    end
end

% generate output text
if foredit
    str = ['% ' str{1} newline];
else
    str = '';
end

if foredit
    for ii = 1:numel(sec)
        if ismember(sec(ii).idx,[2 3])
            str = [str '% ' newline '% ### ' sec(ii).title newline sprintf('%% %s\n',sec(ii).text{:})];
        elseif isempty(sec(ii).text)
            str = [str '% ' newline '% ### ' sec(ii).title newline '% ' newline '% ' newline];
        else
            str = [str '% ' newline '% ### ' sec(ii).title newline '% ' newline sprintf('%% %s\n',sec(ii).text{:})];
        end
    end
else
    for ii = 1:numel(sec)
        if ismember(sec(ii).idx,[2 3])
            str = [str newline '### ' sec(ii).title newline sprintf('%s\n',sec(ii).text{:})];
        elseif isempty(sec(ii).text)
            str = [str newline '### ' sec(ii).title newline newline newline];
        else
            str = [str newline '### ' sec(ii).title newline newline sprintf('%s\n',sec(ii).text{:})];
        end
    end
end
% copy to clipboard
if foredit
    clipboard('copy',str);
end
end