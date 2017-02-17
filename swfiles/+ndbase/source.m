function [dataStr, info] = source(dataSource, type)
% reads text data from multiple sources
%
% [dataStr, {info}] = SOURCE(dataSource, {type})
%
% The function can read text data from the following sources:
%   - local files, where the file location is specified in dataSource
%   - download an url content, when dataSource begins with 'http'
%   - takes the input and pipes it to the output
%
% Input:
%
% dataSource    String, can be file name, url (starting with 'http') or
%               arbitrary string.
% {type}        Forces to treat data source as a certain type, optional.
%               Possible values:
%                   'auto'      Determine automatically, default value.
%                   'file'      Local file.
%                   'url'       Url adress.
%                   'string'    String.
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
    case 'string'
        type = 4;
    otherwise
        error('ndbase:source:WrongInput','The given data source type is invalid!')
end

if type==2 || (type==1 && exist(dataSource,'file')==2)
    % get the name of the file
    fid = fopen(dataSource);
    source = fopen(fid);
    fclose(fid);
    % split the string at new lines
    dataStr = fileread(dataSource);
    isfile = true;
    
elseif type==3 || (type==1 && numel(dataSource)>4 && strcmp('http',dataSource(1:4)) && numel(dataSource)<400)
    % try to load it from the web
    try
        dataStr = urlread(dataSource);
    catch
        error('ndbase:source:WrongInput','The requested remote data is not available!')
    end
    
    source = dataSource;
    isfile = false;
elseif type==4 || type==1
    % just use the input as a string
    dataStr = dataSource;
    
    source = '';
    isfile = false;
end

dataStr = char(dataStr(:))';

if nargout>1
    info.source = source;
    info.isfile = isfile;
end

end