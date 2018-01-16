function export(varargin)
% exports swplot figure into raster/vector image
% 
% ### Syntax
% 
% `swplot.export(Name,Value)`
% 
% ### Description
% 
% `swplot.export(Name,Value)` exports an swplot figure into a raster/vector
% image. The function will remove the tooltip before exporting to raster
% image. For vector graphics, also the legend will be removed as it causes
% a bug in the Matlab [matlab.print] command. Also, vector graphics export
% does not support transparency, thus all transparency will be removed from
% the figure. Be careful, the vector image filesize can be quite large if
% there are many object on the figure. To reduce the file size, try
% reducing the $n_{patch}$ and $n_{mesh}$ values to reduce the number of
% faces per object. The function uses the Matlab built-in [matlab.print]
% command after preparing the figure. All figure property restored after
% export.  
% 
% ### Name-Value Pair Arguments
% 
% `'figure'`
% : Handle of the swplot figure. Default value is the active figure.
% 
% `'filename'`
% : String, name of the image file. Image type will be determined
%   based on the extension. Supported graphics formats:
%   * `png`    Raster image.
%   * `eps`    Vector image.
%
%   If no filename provided, the function returns without printing.
% 
% `'res'`
% : Resolution for raster images in dpi, default value is 300. Set
%   it to 0, to save the image with the screen resolution.
% 
% ### See Also
% 
% [matlab.print]
%

inpForm.fname  = {'figure' 'res' 'filename'};
inpForm.defval = {[]       300   ''        };
inpForm.size   = {[1 1]    [1 1] [1 -1]    };
inpForm.soft   = {true     false true      };

param = sw_readparam(inpForm, varargin{:});

if isempty(param.figure)
    hFigure  = swplot.activefigure('plot');
else
    hFigure = param.figure;
end

if isempty(param.filename)
    return
end

if isempty(param.filename)
    % get file name from UI
    [fName, pathName, ~] = uiputfile({'*.png' 'PNG (raster image)'; '*.eps' 'EPS (vector graphics)'},'Export image as');
    if isnumeric(fName) && fName==0
        warning('export:NoFile','File name is missing!')
        return
    end
    fName = [pathName fName];
    
else
    fName = param.filename;
end

[~,~,ext] = fileparts(fName);

% switch tooltip off
tStatus = swplot.tooltip();
lStatus = swplot.legend;

switch ext
    case '.png'
        swplot.tooltip('off')
        print(hFigure,fName,'-dpng',['-r' num2str(round(param.res))],'-noui');
    case '.eps'
        swplot.tooltip('off')
        swplot.legend('off')
        % find the text/line objects on the figure and hide them
        hPatch = findobj(hFigure,'type','patch');
        % keep alpha mapping
        alpha0 = get(hPatch,'FaceAlpha');
        set(hPatch,'FaceAlpha',1)
        % export vector graphics
        print(hFigure,fName,'-depsc','-noui');
        % restore transparency
        set(hPatch,{'FaceAlpha'},alpha0);
    otherwise
        error('export:WrongExtension','The given file type is unsupported!')
end

swplot.tooltip(tStatus)
swplot.legend(lStatus)

end