function sw_help()
% provide help strings in the compiled application
%
% CHANGE THIS TO HELP()
%

% get the location of the deployed application
if isdeployed && ismac
    NameOfDeployedApp = 'MyDeployedApplication'; % do not include the '.app' extension
    [~, result] = system(['top -n100 -l1 | grep ' NameOfDeployedApp ' | awk ''{print $1}''']);
    result=strtrim(result);
    [status, result] = system(['ps xuwww -p ' result ' | tail -n1 | awk ''{print $NF}''']);
    if status==0
        diridx=strfind(result,[NameOfDeployedApp '.app']);
        realpwd=result(1:diridx-2);
    else
        msgbox({'realpwd not set:',result})
    end
else
    realpwd = pwd;
end

end