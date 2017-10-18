function [strcell,nc,str]=XMLtostring(data,strcell,nc)
if(nargin<2), strcell=cell(1,100000); nc=[0 0]; end
nt=nc(2); nc=nc(1);
nl=length(data);
for i=1:nl
    if(isfieldempty(data(i),'name'))
         nc=nc+1; strcell{nc}=repmat('\t',1,nt);
        if(isfieldempty(data(i),'property'))
            nc=nc+1; strcell{nc}=['<' data(i).name ];
        else
           if(isfieldempty(data(i),'string'))
                nc=nc+1; strcell{nc}=['<' data(i).name '>'];
           else
                nc=nc+1; strcell{nc}=['<' data(i).name '>\n'];
           end
        end
    end
    if(isfieldempty(data(i),'property'))
                for k=1:length(data(i).property)
                    nc=nc+1; strcell{nc}=[' ' data(i).property(k).name '="' data(i).property(k).value '"'];
                end
                if((~isfieldempty(data(i),'node'))&&(~isfieldempty(data(i),'string')))
                    nc=nc+1; strcell{nc}='/>';
                else
                    nc=nc+1; strcell{nc}='>';
                end
                if(~isfieldempty(data(i),'string'))
                    nc=nc+1; strcell{nc}='\n';
                end     
    end
    if(isfieldempty(data(i),'node'))
            [strcell,nc]=XMLtostring(data(i).node,strcell,[nc nt+1]);
    end
    if(isfieldempty(data(i), 'string'))
            nc=nc+1; strcell{nc}=data(i).string;
    end
    if(isfieldempty(data(i),'name'))
        if(isfieldempty(data(i),'node')||isfieldempty(data(i),'string'))
            if(~isfieldempty(data(i),'string'))
                nc=nc+1; strcell{nc}=repmat('\t',1,nt);
            end
            nc=nc+1; strcell{nc}=['</' data(i).name '>\n'];
        end
    end
    if(~(isfieldempty(data(i),'node')||isfieldempty(data(i),'string')||isfieldempty(data(i),'property')))
        nc=nc+1; strcell{nc}=repmat('\t',1,nt);
        nc=nc+1; strcell{nc}=['</' data(i).name '>\n'];
    end
end
str=[strcell{1:nc}];
if(isempty(str)), str=''; end

function a=isfieldempty(data,str)
a=isfield(data,str)&&(~isempty(data.(str)));

