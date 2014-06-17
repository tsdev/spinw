function rootdir = sw_rootdir()
% rootdir = SW_ROOTDIR() gives the location of the spinW program.
%
% See also SW.
%

rootdir = mfilename('fullpath');
idx     = strfind(rootdir,filesep);
rootdir = rootdir(1:idx(end-1));

end