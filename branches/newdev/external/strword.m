function out = strword(str, idx, last)
% extract words separated by whitespace from string
%
% out = strword(str, idx, {last})
%
% Input:
%
% str       String input.
% idx       indexes of the words to be extracted.
% last      If true, the last word is given if idx contains an element
%           larger than the word count, if flase an empty string.
%           Default is false.
%
% Output:
%
% out       Cell contains the extracted words in the order idx is give.
%
% Example:
%
% strword(' one two three four',[4 1])
% the output will be: {'four' 'one'}.
%

if nargin < 3
    last = false;
end

str = strsplit(str);
wordIdx = find(cellfun(@numel,str));

out = cell(1,numel(idx));

for ii = 1:numel(idx)
    
    if wordIdx(end) < idx(ii)
        if last
            out{ii} = str{wordIdx(end)};
        else
            out{ii} = '';
        end
    else
        idxSel = idx(ii);
        out{ii} = str{wordIdx(idxSel)};
    end
    
end

end
