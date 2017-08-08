function [dataStr, info] = source(dataSource, type)
% reads text data from multiple sources
%
% [dataStr, {info}] = SOURCE(dataSource, {type})
%
% The function can read text data from the following sources:
%   - local files, where the file location is specified in dataSource
%   - download an url content, when dataSource begins with 'http'
%
% Input:
%
% dataSource    String, can be file name or url (starting with 'http').
% {type}        Forces to treat data source as a certain type, optional.
%               Possible values:
%                   'auto'      Determine automatically whether file or
%                               url, (default).
%                   'file'      Local file.
%                   'url'       Url adress.
%
% Output:
%
% dataStr       String in a char type row vector.
% {info}        Additional information stored in a struct with fields:
%   source      Data source (empty if string data was given).
%   isfile      Logical variable, true if the source is a local file.
%

if nargin < 2
    type = 'auto';
end

switch type
    case 'auto'
        type = 1;
    case 'file'
        type = 2;
    case 'url'
        type = 3;
    otherwise
        error('ndbase:source:WrongInput','The given data source type is invalid!')
end

if type == 1
    if  exist(dataSource,'file')==2
        % file
        type = 2;
    elseif numel(dataSource)>4 && strcmp('http',dataSource(1:4)) && numel(dataSource)<2000
        % url
        type = 3;
    else
        sourceError(dataSource,type)
    end
end

switch type
    case 2
        % get the name of the file
        fid = fopen(dataSource);
        if fid == -1
            sourceError(dataSource,type)
        end
        source = fopen(fid);
        fclose(fid);
        % split the string at new lines
        dataStr = fileread(dataSource);
        isfile = true;
        
    case 3
        % try to load it from the web
        try
            dataStr = urlread(dataSource);
        catch
            sourceError(dataSource,type)
        end
        
        source = dataSource;
        isfile = false;
end

dataStr = char(dataStr(:))';

if nargout>1
    info.source = source;
    info.isfile = isfile;
end

end

function sourceError(dataSource,type)
% generate error message

if numel(dataSource)<20
    dString = dataSource;
else
    dString = [dataSource(1:8) ' ... ' dataSource(end+(-6:0))];
end

switch type
    case 1
        error('ndbase:source:SourceNotFound','The given data source ''%s'' is neither a valid file nor an existing url!',dString);
    case 2
        error('ndbase:source:SourceNotFound','The given data source ''%s'' is not a valid file!',dString);
    case 3
        error('ndbase:source:SourceNotFound','The given data source ''%s'' is not an existing url!',dString);
end
end