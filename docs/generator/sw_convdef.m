function text0 = sw_convdef(text0,idx,foredit)
% convert list into MD definitions

if nargin<2
    idx = 2;
end

if nargin<3
    foredit = true;
end

% convert arguments and options into definitions
%text0 = sec([sec.idx] == ii).text;
% lines of definition
defIdx   = cellfun(@(C)~isempty(C) && C(1)~=' ' ,text0);
%spaceIdx = cellfun(@(C)find(C==' ',1),text0(defIdx));
[spaceIdx, spaceEnd] = regexp(text0(defIdx),'\s+','once');
spaceEnd = round(mean(cell2mat(spaceEnd)));
text0    = [text0;cell(1,numel(text0))];
fdefIdx  = find(defIdx);
for jj = 1:numel(fdefIdx)
    text0{2,fdefIdx(jj)} = [':' text0{1,fdefIdx(jj)}(spaceEnd:end)];
    if idx == 3
        textTemp = [ '`''' text0{1,fdefIdx(jj)}(1:(spaceIdx{jj}-1)) '''`'];
    else
        textTemp = [ '`' text0{1,fdefIdx(jj)}(1:(spaceIdx{jj}-1)) '`'];
    end
    
    if foredit
        text0{1,fdefIdx(jj)} = [newline '% ' textTemp];
    else
        text0{1,fdefIdx(jj)} = [newline textTemp];
    end
end
text0(1,~defIdx) = cellfun(@(C)[' ' C(spaceEnd:end)],text0(1,~defIdx),'UniformOutput',false);
text0 = text0(cellfun(@(C)~isempty(C),text0(:)));

end