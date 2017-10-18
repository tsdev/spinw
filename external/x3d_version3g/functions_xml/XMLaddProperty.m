function data=XMLaddProperty(name,value,data)
n1=length(data.node);

if(isfieldempty(data.node(n1),'node'))
    n2=length(data.node(n1).node);
    if(isfieldempty(data.node(n1).node(n2),'node'))
        n3=length(data.node(n1).node(n2).node);
        if(isfieldempty(data.node(n1).node(n2).node(n3),'node'))
            n4=length(data.node(n1).node(n2).node(n3).node);
            if(isfieldempty(data.node(n1).node(n2).node(n3).node(n4),'node'))
                n5=length(data.node(n1).node(n2).node(n3).node(n4).node);
                if(isfieldempty(data.node(n1).node(n2).node(n3).node(n4).node(n5),'node'))
                    n6=length(data.node(n1).node(n2).node(n3).node(n4).node(n5).node);
                    if(isfieldempty(data.node(n1).node(n2).node(n3).node(n4).node(n5).node(n6),'node'))
                        n7=length(data.node(n1).node(n2).node(n3).node(n4).node(n5).node(n6).node);
                        if(isfieldempty(data.node(n1).node(n2).node(n3).node(n4).node(n5).node(n6).node(n7),'node'))
                            n8=length(data.node(n1).node(n2).node(n3).node(n4).node(n5).node(n6).node(n7).node);
                            if(isfieldempty(data.node(n1).node(n2).node(n3).node(n4).node(n5).node(n6).node(n7).node(n8),'node'))
                                n9=length(data.node(n1).node(n2).node(n3).node(n4).node(n5).node(n6).node(n7).node(n8).node);
                                if(isfieldempty(data.node(n1).node(n2).node(n3).node(n4).node(n5).node(n6).node(n7).node(n8).node(n9),'node'))
                                    n10=length(data.node(n1).node(n2).node(n3).node(n4).node(n5).node(n6).node(n7).node(n8).node(n9).node);
                                    if(isfieldempty(data.node(n1).node(n2).node(n3).node(n4).node(n5).node(n6).node(n7).node(n8).node(n9).node(n10),'node'))
                                        disp('reached maximum recursion limit');
                                    else
                                        if(~isfieldempty(data.node(n1).node(n2).node(n3).node(n4).node(n5).node(n6).node(n7).node(n8).node(n9).node(n10),'property')), np=1; else np=length(data.node(n1).node(n2).node(n3).node(n4).node(n5).node(n6).node(n7).node(n8).node(n9).node(n10).property)+1; end
                                        data.node(n1).node(n2).node(n3).node(n4).node(n5).node(n6).node(n7).node(n8).node(n9).node(n10).property(np).name=name; data.node(n1).node(n2).node(n3).node(n4).node(n5).node(n6).node(n7).node(n8).node(n9).node(n10).property(np).value=value;
                                    end
                                else
                                    if(~isfieldempty(data.node(n1).node(n2).node(n3).node(n4).node(n5).node(n6).node(n7).node(n8).node(n9),'property')), np=1; else np=length(data.node(n1).node(n2).node(n3).node(n4).node(n5).node(n6).node(n7).node(n8).node(n9).property)+1; end
                                    data.node(n1).node(n2).node(n3).node(n4).node(n5).node(n6).node(n7).node(n8).node(n9).property(np).name=name; data.node(n1).node(n2).node(n3).node(n4).node(n5).node(n6).node(n7).node(n8).node(n9).property(np).value=value;
                                end
                            else
                                if(~isfieldempty(data.node(n1).node(n2).node(n3).node(n4).node(n5).node(n6).node(n7).node(n8),'property')), np=1; else np=length(data.node(n1).node(n2).node(n3).node(n4).node(n5).node(n6).node(n7).node(n8).property)+1; end
                                data.node(n1).node(n2).node(n3).node(n4).node(n5).node(n6).node(n7).node(n8).property(np).name=name; data.node(n1).node(n2).node(n3).node(n4).node(n5).node(n6).node(n7).node(n8).property(np).value=value;
                            end
                        else
                            if(~isfieldempty(data.node(n1).node(n2).node(n3).node(n4).node(n5).node(n6).node(n7),'property')), np=1; else np=length(data.node(n1).node(n2).node(n3).node(n4).node(n5).node(n6).node(n7).property)+1; end
                            data.node(n1).node(n2).node(n3).node(n4).node(n5).node(n6).node(n7).property(np).name=name; data.node(n1).node(n2).node(n3).node(n4).node(n5).node(n6).node(n7).property(np).value=value;
                        end
                    else
                        if(~isfieldempty(data.node(n1).node(n2).node(n3).node(n4).node(n5).node(n6),'property')), np=1; else np=length(data.node(n1).node(n2).node(n3).node(n4).node(n5).node(n6).property)+1; end
                        data.node(n1).node(n2).node(n3).node(n4).node(n5).node(n6).property(np).name=name; data.node(n1).node(n2).node(n3).node(n4).node(n5).node(n6).property(np).value=value;
                    end
                else
                    if(~isfieldempty(data.node(n1).node(n2).node(n3).node(n4).node(n5),'property')), np=1; else np=length(data.node(n1).node(n2).node(n3).node(n4).node(n5).property)+1; end
                    data.node(n1).node(n2).node(n3).node(n4).node(n5).property(np).name=name; data.node(n1).node(n2).node(n3).node(n4).node(n5).property(np).value=value;
                end
            else
                if(~isfieldempty(data.node(n1).node(n2).node(n3).node(n4),'property')), np=1; else np=length(data.node(n1).node(n2).node(n3).node(n4).property)+1; end
                data.node(n1).node(n2).node(n3).node(n4).property(np).name=name; data.node(n1).node(n2).node(n3).node(n4).property(np).value=value;
            end
        else
            if(~isfieldempty(data.node(n1).node(n2).node(n3),'property')), np=1; else np=length(data.node(n1).node(n2).node(n3).property)+1; end
            data.node(n1).node(n2).node(n3).property(np).name=name; data.node(n1).node(n2).node(n3).property(np).value=value;
        end
    else
        if(~isfieldempty(data.node(n1).node(n2),'property')), np=1; else np=length(data.node(n1).node(n2).property)+1; end
        data.node(n1).node(n2).property(np).name=name; data.node(n1).node(n2).property(np).value=value;
    end
else
    if(~isfieldempty(data.node(n1),'property')), np=1; else np=length(data.node(n1).property)+1; end
    data.node(n1).property(np).name=name; data.node(n1).property(np).value=value;
end

function a=isfieldempty(data,str)
a=isfield(data,str)&&(~isempty(data.(str)));