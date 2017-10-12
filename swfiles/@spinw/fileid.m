function fidOut = fileid(obj,fid)
% determines file object for text output
% 
% ### Syntax
% 
% `fileid(obj,fid)`
% `fidOut = fileid(obj)`
% 
% ### Description
% 
% `fileid(obj,fid)` determines the text output of all [spinw] class
% methods. Default is 1, where all output is printed onto the MATLAB
% Command Window.
%  
% `fidOut = fileid(obj)` outputs the stored fileID value.
% 
% ### See Also
% 
% [spinw] \| [spinw.disp]
%

if nargin > 1
    obj.fid = fid;
    return
end

fidOut = obj.fid;

end