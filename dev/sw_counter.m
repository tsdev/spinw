function sw_counter(reset, text)
% prints the number of function calls
% 
% ### Syntax
% 
% `sw_counter`
%
% `sw_counter(reset, text)`
% 
% ### Description
% 
% `sw_counter(reset, text)` prints the number of calls to this function.
%
% ### Input Arguments
%
% `reset`
% : If `true` the counter is reset to zero.
%
% `text`
% : Sets the text to print before the number, default value is `'The number
%   of function calls:'`.
%

persistent cnt

if nargin == 0
    reset = 0;
end

if reset
    cnt = 0;
else
    if isempty(cnt)
        cnt = 0;
    end
    cnt = cnt + 1;
end

if nargin < 2
    text = 'The number of function calls:';
end

nT = numel(text);
bb = repmat('\b',[1 nT+7]);

if cnt == 0
    fprintf(['\n' text '%6d\n'],cnt);
else
    fprintf([bb text '%6d\n'],cnt);
end

end