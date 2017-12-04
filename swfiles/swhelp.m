function varargout = swhelp(funName0)
% outputs the SpinW help
%
% ### Syntax
%
% `swhelp(funName)`
%
% ### Description
%
% `swhelp(funName)` shows the help on the given function name,
% method, property name. Works the same way as the Matlab built-in
% [matlab.help] command.
%
% ### See Also
%
% [swdoc], [matlab.help]
%

helpStr = help(funName0);
% remove the  --- method --- part
helpStr = regexprep(helpStr,'---[^\n]+?---\n\s*?\n','');
newLine = char(10); %#ok<CHARTEN>

% if there is a ### Methods section convert all into links
mIdx1 = regexp(helpStr,'### Methods','once');
if ~isempty(mIdx1)
    mIdx2 = regexp(helpStr,' ### ');
    mIdx2 = mIdx2(find(mIdx2>mIdx1,1));
    
    % replace links
    %helpStr = [helpStr(1:mIdx1-1) regexprep(helpStr(mIdx1:mIdx2),['(' funName0 '\.\w+?)(\s+)'],'\[$1\] ${strtrim(iindex(strsplit(help($1),char(10)),''{}'',1))}$2') helpStr(mIdx2+1:end)];
    helpStr = [helpStr(1:mIdx1-1) regexprep(helpStr(mIdx1:mIdx2),['[ ]+(' funName0 '\.\w+?)\n'],'  \* \[$1\] ${strtrim(iindex(strsplit(help($1),char(10)),''{}'',1))}\n') helpStr(mIdx2+1:end)];
end

% remove reference page thingy
helpStr = regexprep(helpStr,'Reference page in Doc Center\n[\s\w]+?\n','');

helpStr = [repmat('_',1,75) newline '`' funName0 '` ' helpStr];
    
% convert from SpinW MarkDown to output that can be printed in the Command Window
helpStr = sw_markdown(helpStr);

% cut trailing space
helpStr = regexprep(helpStr,'\s*$','');

helpStr = [helpStr newLine repmat('_',1,75)];


if nargout == 0
    % print the help
    disp(helpStr);
else
    varargout{1} = helpStr;
end

end