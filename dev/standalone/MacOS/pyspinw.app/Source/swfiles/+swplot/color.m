function [RGB, nameOut] = color(cName, index)
% generates RGB code from color name string
%
% RGB = SWPLOT.COLOR(cName,{index})
%
% Input:
%
% cName     String that contains the name of the color, either a single
%           character (see <a href="matlab: doc ColorSpec">ColorSpec</a>) or use any HTML color name,
%           (see http://www.w3schools.com/html/html_colornames.asp).
%           For multiple colors, use a cell containing the strings. The
%           name of the colors are stored in the <a href="matlab: edit color.dat">color.dat</a> file.
% index     If true, the index of the color in the color.dat file is read.
%           Index 1 corresponds to the 9th entry (the first 8 entry has
%           already names in Matlab). Default is false.
%
% Output:
%
% RGB       RGB color code, dimensions are [3 nColor], where
%           every value is between 0 and 255.
%
% Example:
%
%   RGB = SWPLOT.COLOR('LightGray')
%
%   the output RGB will be [211; 211; 211].
%

if nargin == 0
    help swplot.color
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