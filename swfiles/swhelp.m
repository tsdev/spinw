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

% convert []() links into anchor tag
if feature('hotlinks')
    helpStr = regexprep(helpStr,'\[([^\n]+?)\]\((.+?)\)','<a href="$2">$1</a>');
else
    % no links
    helpStr = regexprep(helpStr,'\[([^\n]+?)\]\((.+?)\)','$1');
end

% remove abbreviation list
helpStr = regexprep(helpStr,'\n\s*?\*\[[.a-zA-Z0-9\s]+?\]:.+?\n','\n');
% remove code formatting
helpStr = regexprep(helpStr,'>>>','');
helpStr = regexprep(helpStr,'>>','');
% remove snapnow code
helpStr = regexprep(helpStr,'\n\s*snapnow\n','\n');
helpStr = regexprep(helpStr,'\\\\(\w+)','${symbol($1,2)}');

regexp0 = '\[\!\[([^\!\[\]]+?)\]\((\S+?)\)\{(.+?)\}\]\((\S+?)\)';
helpStr = regexprep(helpStr,regexp0,'[image "$2", alt="$1", caption="$3"]');
% image with caption only
regexp0 = '\!\[([^\!\[\]]+?)\]\((\S+?)\)\{(.+?)\}';
helpStr = regexprep(helpStr,regexp0,'[image "$2", alt="$1", caption="$3"]');
% image without caption
regexp0 = '\!\[([^\!\[\]]+?)\]\((\S+?)\)';
helpStr = regexprep(helpStr,regexp0,'[image "$2", alt="$1"]');

% link to open .dat files in editor
if feature('hotlinks')
    helpStr = regexprep(helpStr,'\[(\w+?)\.dat\]','<a href="matlab:edit([sw_rootdir,''dat_files'',filesep,''$1.dat''])">$1.dat</a>');
    % create links
    helpStr = regexprep(helpStr,'\[matlab\.(\w+?)\]','<a href="matlab:help(''$1'')">$1</a>');
    helpStr = regexprep(helpStr,'\[(\w+?)\.(\w+?)\]','<a href="matlab:swhelp(''$1.$2'')">$1.$2</a>');
    helpStr = regexprep(helpStr,'\[(\w+?)\]','<a href="matlab:swhelp(''$1'')">$1</a>');
    
else
    helpStr = regexprep(helpStr,'\[(\w+?)\.dat\]','$1.dat');
    
    helpStr = regexprep(helpStr,'\[matlab\.(\w+?)\]','$1');
    helpStr = regexprep(helpStr,'\[(\w+?)\.(\w+?)\]','$1.$2');
    helpStr = regexprep(helpStr,'\[(\w+?)\]','$1');

end

% remove the note marking
helpStr = regexprep(helpStr,'{{warning (.+?)}}',[newline repmat([symbol('skull') ' '],1,37) newline '$1' newline repmat([symbol('skull') ' '],1,37)]);
helpStr = regexprep(helpStr,'{{note (.+?)}}',[newline repmat('~',1,75) newline '$1' newline repmat('~',1,75)]);

% cut trailing space
helpStr = regexprep(helpStr,'\s*$','');

helpStr = [repmat('_',1,75) newline '`' funName0 '` ' helpStr newLine repmat('_',1,75)];

if feature('hotlinks')
    helpStr = regexprep(helpStr,'`(\s)','</strong>$1');
    helpStr = regexprep(helpStr,'(\s)`','$1<strong>');
    helpStr = regexprep(helpStr,'`([\n\,\.\)])','</strong>$1');
    helpStr = regexprep(helpStr,'([\n\(])`','$1<strong>');
    helpStr = regexprep(helpStr,'\\mathbf\{(\w+?)\}','<strong>$1</strong>');
else
    helpStr = regexprep(helpStr,'`','');
    helpStr = regexprep(helpStr,'\\mathbf\{(\w+?)\}','$1');
end
helpStr = regexprep(helpStr,'`','');
helpStr = regexprep(helpStr,'\$','');
helpStr = regexprep(helpStr,'\\begin\{align\}','             ');
helpStr = regexprep(helpStr,'\\end\{align\}','           ');
helpStr = regexprep(helpStr,' \\\|',',');
helpStr = regexprep(helpStr,'\s*\\times\s*',symbol('cross'));
helpStr = regexprep(helpStr,'\*\*','');
helpStr = regexprep(helpStr,'_\{(\w+?)\}','_$1');
helpStr = regexprep(helpStr,'\^\{(\w+?)\}','^$1');
helpStr = regexprep(helpStr,'\&=','=');
helpStr = regexprep(helpStr,'\\\\','');
helpStr = regexprep(helpStr,'\^\{-1\}',[symbol('^-') symbol('^1')]);
helpStr = regexprep(helpStr,'\s*\\cdot\s*',symbol('cross'));

helpStr = regexprep(helpStr,'\s*\\rangle',symbol('ket'));
helpStr = regexprep(helpStr,'\\langle\s*',symbol('bra'));

helpStr = regexprep(helpStr,'\\(\w+)','${symbol($1,2)}');
helpStr = regexprep(helpStr,'\^2',symbol('^2'));
helpStr = regexprep(helpStr,'\^3',symbol('^3'));
%helpStr = regexprep(helpStr,'\n\s*?\n\s*?\n',[newLine newLine]);


if nargout == 0
    % print the help
    disp(helpStr);
else
    varargout{1} = helpStr;
end

end