function [data,loc]=XMLaddNode(name,data,loc)
if(nargin<2)
    data=struct(); n1=0;
else
    if(isfieldempty(data,'node')), n1=length(data.node); else n1=0; end
end
if(nargin<3)
    loc=0;
end
if(loc==0);
    data.node(n1+1).name=name;
else
    if(isfieldempty(data.node(n1),'node')), n2=length(data.node(n1).node); else n2=0; end
    if(loc==1)
        data.node(n1).node(n2+1).name=name;
    else
        if(isfieldempty(data.node(n1).node(n2),'node')), n3=length(data.node(n1).node(n2).node); else n3=0; end
        if(loc==2)
            data.node(n1).node(n2).node(n3+1).name=name;
        else
            if(isfieldempty(data.node(n1).node(n2).node(n3),'node')), n4=length(data.node(n1).node(n2).node(n3).node); else n4=0; end
            if(loc==3)
                data.node(n1).node(n2).node(n3).node(n4+1).name=name;
            else
                if(isfieldempty(data.node(n1).node(n2).node(n3).node(n4),'node')), n5=length(data.node(n1).node(n2).node(n3).node(n4).node); else n5=0; end
                if(loc==4)
                    data.node(n1).node(n2).node(n3).node(n4).node(n5+1).name=name;
                else
                    if(isfieldempty(data.node(n1).node(n2).node(n3).node(n4).node(n5),'node')), n6=length(data.node(n1).node(n2).node(n3).node(n4).node(n5).node); else n6=0; end
                    if(loc==5)
                        data.node(n1).node(n2).node(n3).node(n4).node(n5).node(n6+1).name=name;
                    else
                        if(isfieldempty(data.node(n1).node(n2).node(n3).node(n4).node(n5).node(n6),'node')), n7=length(data.node(n1).node(n2).node(n3).node(n4).node(n5).node(n6).node); else n7=0; end
                        if(loc==6)
                            data.node(n1).node(n2).node(n3).node(n4).node(n5).node(n6).node(n7+1).name=name;
                        else
                            if(isfieldempty(data.node(n1).node(n2).node(n3).node(n4).node(n5).node(n6).node(n7),'node')), n8=length(data.node(n1).node(n2).node(n3).node(n4).node(n5).node(n6).node(n7).node); else n8=0; end
                            if(loc==7)
                                data.node(n1).node(n2).node(n3).node(n4).node(n5).node(n6).node(n7).node(n8+1).name=name;
                            else
                                if(isfieldempty(data.node(n1).node(n2).node(n3).node(n4).node(n5).node(n6).node(n7).node(n8),'node')), n9=length(data.node(n1).node(n2).node(n3).node(n4).node(n5).node(n6).node(n7).node(n8).node); else n9=0; end
                                if(loc==8)
                                    data.node(n1).node(n2).node(n3).node(n4).node(n5).node(n6).node(n7).node(n8).node(n9+1).name=name;
                                else
                                        if(isfieldempty(data.node(n1).node(n2).node(n3).node(n4).node(n5).node(n6).node(n7).node(n8).node(n9),'node')), n10=length(data.node(n1).node(n2).node(n3).node(n4).node(n5).node(n6).node(n7).node(n8).node(n9).node); else n10=0; end
                                        if(loc==9)
                                            data.node(n1).node(n2).node(n3).node(n4).node(n5).node(n6).node(n7).node(n8).node(n9).node(n10+1).name=name;
                                        else
                                            disp('reached maximum recursion limit');
                                        end
                                end

                            end
                        end
                    end
                end
            end
        end
    end
end

function a=isfieldempty(data,str)
a=isfield(data,str)&&(~isempty(data.(str)));

