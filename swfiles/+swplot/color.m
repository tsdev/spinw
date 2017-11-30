function [RGB, nameOut] = color(cName, index)
% generates RGB code from color name
% 
% ### Syntax
% 
% `RGB = swplot.color(cName)`
%
% `RGB = swplot.color(cName,index)`
% 
% ### Description
% 
% `RGB = swplot.color(cName)` reads the color RGB values from the
% [color.dat] file corresponding to the given color name `cName`. The
% color name can be either a single character (see [matlab.colorspec]) or
% any [HTML color name](https://www.w3schools.com/colors/colors_names.asp)
% that is stored in the [color.dat] file.
% 
% `RGB = swplot.color(cName,index)` if `index` is true, RGB code
% corresponding to the `cName` color index is read.
%
% ### Examples
% 
% Read the RGB code corresponding to light gray:
% ```
% >>RGB = swplot.color('LightGray')>>
% ```
% 
% ### Input Arguments
% 
% `cName`
% : String of a color name. For multiple colors, use a cell of strings.
% 
% `index`
% : If `true`, instead of the color name, `cName` means the index of the
%   color in the [color.dat] file. index 1 corresponds to the 9th entry
%   (the first 8 entry are standard Matlab color names), default value is
%   `false`.
% 
% ### Output Arguments
% 
% `RGB`
% : RGB color codes in a matrix with dimensions of $[3\times n_{color}]$, where
%   every value is an integer between 0 and 255.
%

if nargin == 0
    swhelp swplot.color
    return
end

if nargin < 2
    index = false;
end

if index
    % Open the color definition file.
    colPath = [sw_rootdir 'dat_files' filesep 'color.dat'];
    fid = fopen(colPath);
    if fid == -1
        error('color:FileNotFound',['Color definition file not found: '...
            regexprep(colPath,'\' , '\\\') '!']);
    end
    % read color names and RGB values from the file
    cDat = textscan(fid,'%s %d %d %d');
    cNameList = cDat{1};
    RGBList = [cDat{2:4}];
    
    fclose(fid);
    
    RGB = double(RGBList(cName+8,:)');
    if nargout > 1
        nameOut = cNameList(cName+8);
    end
    
    return
end

if isnumeric(cName)
    if ~any(size(cName)-[1 3])
        RGB = cName';
        return
    end
    
    if size(cName,1) == 3
        RGB = cName;
        return
    end
    
    RGB = [];
end

if ischar(cName)
    cName = {cName};
end

if iscell(cName)
    % Open the color definition file.
    colPath = [sw_rootdir 'dat_files' filesep 'color.dat'];
    fid = fopen(colPath);
    if fid == -1
        error('color:FileNotFound',['Color definition file not found: '...
            regexprep(colPath,'\' , '\\\') '!']);
    end
    % read color names and RGB values from the file
    cDat = textscan(fid,'%s %d %d %d');
    cNameList = cDat{1};
    RGBList = [cDat{2:4}];
    
    fclose(fid);
    
    % loop though for all color names
    idx = zeros(1,numel(cName));
    for ii = 1:numel(cName)
        % find color index in the list
        idxtmp = find(strcmpi(cNameList,cName{ii}),1);
        if ~isempty(idxtmp)
            idx(ii) = idxtmp;
        end
    end
    
    if any(idx==0)
        RGB = [];
        eIdx = find(idx==0,1);
    else
        RGB = RGBList(idx,:);
    end
else
    cName = {cName};
end

if isempty(RGB)
    error('color:WrongInput','The given color name (''%s'') does not exists!',cName{eIdx})
end

RGB = double(RGB');

end