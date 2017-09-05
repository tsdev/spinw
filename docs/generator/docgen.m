%% generate help

sw_genhelp('path',[sw_rootdir 'swfiles'])
sw_genhelp('path',[sw_rootdir 'swfiles/@spinw'])
sw_genhelp('path',[sw_rootdir 'swfiles/+swplot'])
sw_genhelp('path',[sw_rootdir 'swfiles/+swpref'])
sw_genhelp('path',[sw_rootdir 'swfiles/+swsym'])

%%
sidebar = struct;
sidebar.entries.title   = 'Sidebar';
sidebar.entries.product = 'SpinW';
sidebar.entries.version = '3.0.1';
sidebar.entries.folder(1).title  = '';
sidebar.entries.folder(1).output = 'pdf';
sidebar.entries.folder(1).type   = 'frontmatter';
sidebar.entries.folder(1).folderitems(1).title = '';
sidebar.entries.folder(1).folderitems(1).url   = '/titlepage.html';

YAML.dump(sidebar)