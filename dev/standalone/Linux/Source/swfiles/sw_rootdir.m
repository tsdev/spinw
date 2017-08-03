function rootdir = sw_rootdir()
% gives the path to the SpinW code
%
% rootdir = SW_ROOTDIR()
%
% See also SPINW.
%

rootdir = mfilename('fullpath');
idx     = strfind(rootdir,filesep);
rootdir = rootdir(1:idx(end-1));

end