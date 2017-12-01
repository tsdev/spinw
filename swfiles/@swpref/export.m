function success = export(obj,location)
% saves swpref object into a file
%
% ### Syntax
%
% `success = export(obj,location)`
%
% `success = export(obj)`
%
% ### Description
%
% `success = export(obj,location)` writes the preferences given in `obj` to
% a file location given by `location`. The file is in a basic `.json`
% format.
%
% `success = export(obj)` writes the preferences given in `obj` to
% the users home folder as `swprefs.json`. The file is in a basic `.json`
% format.
%
% ### See Also
%
% [swpref.import]
%

props = obj.props;

if nargin == 1
    location = [userpath filesep 'swprefs.json'];
end


f = fopen(location,'w');
if f < 0
    % error
    error('swpref:export:FileBlocked','Cannot write to file %s.',location);
end
fprintf(f, '{\n');
for i = 1:length(props)
    name = props(i).Name;
    value = obj.(name);
    if ~isstruct(value)
        if isnumeric(value)
            if i == length(props)
                fprintf(f, '\t''%s'': %f\n',name,value);
            else
                fprintf(f, '\t''%s'': %f,\n',name,value);
            end
        elseif ischar(value)
            if i == length(props)
                fprintf(f, '\t''%s'': ''%s''\n',name,value);
            else
                fprintf(f, '\t''%s'': ''%s'',\n',name,value);
            end
        elseif isa(value,'function_handle')
            if i == length(props)
                fprintf(f, '\t''%s'': ''@%s''\n',name,func2str(value));
            else
                fprintf(f, '\t''%s'': ''@%s'',\n',name,func2str(value));
            end
        end
    else
    end
end
fprintf(f, '}\n');
success = fclose(f);

end