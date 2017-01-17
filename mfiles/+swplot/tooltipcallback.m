function tooltipcallback(obj,hit,hFigure,hTransform)
% call tooltip on clicking an object on an swplot figure
%
% SWPLOT.TOOLTIPCALLBACK(obj,hit,hFigure,hTransform)
%
% Input:
%
% obj       Object that was clicked on.
% hit       Hit object, defines the point where the mouse clicked.
% hFigure   Handle of parent swplot figure.
% hTransform Parent hgtransform object if exists.
%
% See also SWPLOT.TOOLTIP.
%


if isappdata(obj,'facenumber')
    % face patch object
    % convert face into object index
    fIdx = swplot.patchfacefcn(obj,hit,[],'all',[1e-7 1e-7],hTransform);
    
    if ~isempty(fIdx)
        fn   = getappdata(obj,'facenumber');
        % object number
        number = fn(fIdx(1));
        
        % graphical object data on swplot figure
        sObject = getappdata(hFigure,'objects');
        
        % find the corresponding object data index
        nIdx = find([sObject(:).number]==number);
        % find text and show tooltip if there is text
        if ~isempty(nIdx)
            % selected data object
            sObject = sObject(nIdx);
            if ismember(sObject.name,{'atom' 'mag'})
                % generate text automatically
                string = swplot.tooltipstring(sObject);
            else
                % use the given text
                string  = sObject.text;
            end
            
            if isempty(string)
                % add some basic information
                if strcmp(sObject.type,'facepatch')
                    sObject.type = 'patch';
                end
                string = [sObject.type ' ' num2str(sObject.number)];
                
                
                % generate position text
                posVal = ~any(isnan(sObject.position),1);
                pos = sObject.position(:,posVal);
                
                
                if ~isempty(pos)
                    posLabel = {'r_1' 'r_2'};
                    posLabel = posLabel(posVal);
                    for ii = 1:size(pos,2)
                        string = [string sprintf(['\n' posLabel{ii} ' = [%5.2f,%5.2f,%5.2f]'],pos(1,ii),pos(2,ii),pos(3,ii))]; %#ok<AGROW>
                    end
                end
            end
            
            swplot.tooltip(string,hFigure);
            
        end
    else
        swplot.tooltip('',hFigure);
    end
elseif isappdata(obj,'vertexnumber')
    % edge object
    swplot.tooltip('',hFigure);
else
    % graphical object data on swplot figure
    sObject = getappdata(hFigure,'objects');
    nIdx = find(obj == [sObject(:).handle]);
    % find text and show tooltip if there is text
    if ~isempty(nIdx)
        swplot.tooltip(sObject(nIdx(1)).text,hFigure)
    end
    
end

end