function [latt_parms,Bmat,proj_out,D,signal,name] = read_struct(dstruct,proj,dproj)
%dstruct is the spinw structure
%proj is the viewing projection matrix
%drproj is the size of the bin in each direction (ignored for the direction
%of the cut.
    fnames=fieldnames(dstruct);
    objnum = find(strcmp(fnames,'obj'));
    swobj = dstruct.(subsref(fnames,substruct('{}',{objnum})));
    latt_parms = abc(swobj);
    M = basisvector(swobj);
    Bmat = inv(M);
    name =formula(swobj).chemform;
    proj_out = proj(:);
    hkls = dstruct.hkl;
    hkls_sz = size(hkls);
    dir_vec=hkls(:,hkls_sz(2))-hkls(:,1);
    D={};
    for idx=1:3
        procjv = proj(idx,:)/norm(proj(idx,:));
       
       if abs(sum(cross(dir_vec,procjv)))< 1e-6
           D{idx} = zeros([1,hkls_sz(2)+1]);
           for idx2=1:hkls_sz(2)
               D{idx}(idx2) = dot(hkls(:,idx2),procjv);
           end
           %move to bin boundaries
           D{idx}(2:hkls_sz(2)) = (D{idx}(1:(hkls_sz(2)-1))+D{idx}(2:(hkls_sz(2))))/2;
           % add a bin boundary on the end
           D{idx}(hkls_sz(2)+1)=2*D{idx}(hkls_sz(2))-D{idx}(hkls_sz(2)-1);
           % add a bin boundary to the beginning
           D{idx}(1) = 2*D{idx}(2)-D{idx}(3);
       else
           hkl_proj = dot(hkls(:,1),procjv);
           D{idx} = [hkl_proj-dproj(idx)/2, hkl_proj+dproj(idx)/2];
       end  
    end
    D{4}=dstruct.Evect;
    signal = dstruct.swConv;
    

