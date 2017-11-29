function success = export(obj,location)
%EXPORT Summary of this function goes here
%   Detailed explanation goes here

props = obj.props;

if nargin == 1
    location = [sw_rootdir filesep 'prefs.json'];
end


f = fopen(location,'w');
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