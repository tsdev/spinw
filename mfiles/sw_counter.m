function sw_counter(reset, text)
% print the number of calls to this functio to the Command Line

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