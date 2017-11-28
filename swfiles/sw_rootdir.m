function rootdir = sw_rootdir
% path to the SpinW folder
% 
% ### Syntax
% 
% `rootdir = sw_rootdir`
% 
% ### Description
% 
% `rootdir = sw_rootdir` returns the parent folder of the `swfiles` folder.
% 
% ### See Also
% 
% [spinw]
%

rootdir = mfilename('fullpath');
idx     = strfind(rootdir,filesep);
rootdir = rootdir(1:idx(end-1));

end