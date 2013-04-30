function sw_initialize(varargin)
% SW_INITIALIZE initializes sw library and creates help files in the Matlab
% documentation.
%
% See also SW.
%

cd(sw_rootdir);
m2html('mfiles','m_files','htmldir','doc')
builddocsearchdb([sw_rootdir 'doc']);

end