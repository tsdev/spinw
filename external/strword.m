function out = strword(str, idx, last, delimiter)
% extract words separated by a delimiter from a string
%
% out = strword(str, idx, {last},{delimiter})
%
% The function extracts multiple words from a string with the given idx
% indices. If a member of idx exceeds the number of words in str, an empty
% string is returned unless last is true. In this case the last word is
% returned.
%
% Input:
%
% str       String input.
% idx       indexes of the words to be extracted.
% last      If true, the last word is returned idx>nWord, if false empty
%           string is returned for idx>wordCount. Default is false.
% delimiter Delimiter that separates words. Special characters can be given
%           with the standard sprintf notations, such as '\t' for tab.
%           Default value is whitespace.
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
if nargin < 4
    delimiter = ' ';
else
    % handle special characters 
    delimiter = sprintf(delimiter);
end

str = strsplit(str,delimiter);
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