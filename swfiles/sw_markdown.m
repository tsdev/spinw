function str = sw_markdown(str,hotlinks)
% converts markdown like text
%
% ### Syntax
%
% `sw_markdown(str,hotlinks)`
%
% ### Description
%
% `sw_markdown(str,hotlinks)` shows the help on the given function name,
% method, property name. Works the same way as the Matlab built-in
% [matlab.help] command.
%
% ### See Also
%
% [swdoc], [matlab.help]
%

if nargin<2
    hotlinks = feature('HotLinks');
end

% convert []() links into anchor tag
if hotlinks
    str = regexprep(str,'\[([^\n^\]]+?)\]\((.+?)\)','<a href="$2">$1</a>');
else
    % no links
    str = regexprep(str,'\[([^\n^\]]+?)\]\((.+?)\)','$1');
end

% remove abbreviation list
str = regexprep(str,'\s*?\*\[[.a-zA-Z0-9\s]+?\]:.+?\n','');
% remove code formatting
str = regexprep(str,'>>>','');
str = regexprep(str,'>>','');
% remove snapnow code
str = regexprep(str,'\n\s*snapnow\n','\n');
str = regexprep(str,'\\\\([a-zA-Z]+)','${symbol($1,2)}');

regexp0 = '\[\!\[([^\!\[\]]+?)\]\((\S+?)\)\{(.+?)\}\]\((\S+?)\)';
str = regexprep(str,regexp0,'[image "$2", alt="$1", caption="$3"]');
% image with caption only
regexp0 = '\!\[([^\!\[\]]+?)\]\((\S+?)\)\{(.+?)\}';
str = regexprep(str,regexp0,'[image "$2", alt="$1", caption="$3"]');
% image without caption
regexp0 = '\!\[([^\!\[\]]+?)\]\((\S+?)\)';
str = regexprep(str,regexp0,'[image "$2", alt="$1"]');

% link to open .dat files in editor
if hotlinks
    str = regexprep(str,['\[([\w^\' filesep ']+?)\.dat\]'],'<a href="matlab:edit([sw_rootdir,''dat_files'',filesep,''$1.dat''])">$1.dat</a>');
    str = regexprep(str,['\[([\w\~\\\/]+?\' filesep '{1,1}[\w]+?[\.]+?[\w]+?)\]'],'<a href="matlab:edit(''$1'')">$1</a>');
    % create links
    str = regexprep(str,'\[matlab\.(\w+?)\]','<a href="matlab:help(''$1'')">$1</a>');
    str = regexprep(str,'\[(\w+?)\.(\w+?)\]','<a href="matlab:swhelp(''$1.$2'')">$1.$2</a>');
    str = regexprep(str,'\[(\w+?)\]','<a href="matlab:swhelp(''$1'')">$1</a>');
    
else
    str = regexprep(str,['\[([\w^\' filesep ']+?)\.dat\]'],'$1.dat');
    str = regexprep(str,['\[([\w\~\\\/]+?\' filesep '{1,1}[\w]+?[\.]+?[\w]+?)\]'],'$1');
    
    str = regexprep(str,'\[matlab\.(\w+?)\]','$1');
    str = regexprep(str,'\[(\w+?)\.(\w+?)\]','$1.$2');
    str = regexprep(str,'\[(\w+?)\]','$1');

end

% remove the note marking
str = regexprep(str,'{{warning (.+?)}}',[newline repmat([symbol('skull') ' '],1,37) newline '$1' newline repmat([symbol('skull') ' '],1,37)]);
str = regexprep(str,'{{note (.+?)}}',[newline repmat('~',1,75) newline '$1' newline repmat('~',1,75)]);

if hotlinks
    str = regexprep(str,'`([\s\,\.\)\:])','</strong>$1');
    str = regexprep(str,'^`','<strong>');
    str = regexprep(str,'`$','</strong>');
    str = regexprep(str,'([\s\(])`','$1<strong>');
    str = regexprep(str,'\\mathbf\{(\w+?)\}','<strong>$1</strong>');
else
    str = regexprep(str,'`','');
    str = regexprep(str,'\\mathbf\{(\w+?)\}','$1');
end
str = regexprep(str,'`','');
str = regexprep(str,'\$','');
str = regexprep(str,'\\begin\{align\}','             ');
str = regexprep(str,'\\end\{align\}','           ');
str = regexprep(str,' \\\|',',');
str = regexprep(str,'\s*\\times\s*',symbol('cross'));
str = regexprep(str,'\*\*','');
str = regexprep(str,'_\{(\w+?)\}','_$1');
str = regexprep(str,'\^\{(\w+?)\}','^$1');
str = regexprep(str,'\&=','=');
str = regexprep(str,'\\\\','');
str = regexprep(str,'\^\{-1\}',[symbol('^-') symbol('^1')]);
str = regexprep(str,'\s*\\cdot\s*',symbol('cross'));

str = regexprep(str,'\s*\\rangle',symbol('ket'));
str = regexprep(str,'\\langle\s*',symbol('bra'));

str = regexprep(str,'\\(\w+)','${symbol($1,2)}');
str = regexprep(str,'([\^\_]{1,1})\-([0-9]{1,1})','${symbol([''\\\\'' $1 ''-\\\\'' $1 $2])}');
str = regexprep(str,'([\^\_]{1,1}[0-9]{1,1})','${symbol($1)}');

%helpStr = regexprep(helpStr,'\n\s*?\n\s*?\n',[newLine newLine]);

end