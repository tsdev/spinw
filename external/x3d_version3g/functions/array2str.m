function str=array2str(A,stre)
    A(isnan(A))=0;
    if(isa(A,'integer'))
        row=repmat('%d ',1,size(A,2));
    else
        row=repmat('%4.5f ',1,size(A,2));
    end

    if(nargin>1)
        row=[row stre];
    end
    str=sprintf(row,A');