function str=interactive_js(number,Tag0)
ii=0;
ii=ii+1; lines{ii}='<![CDATA[\n';
ii=ii+1; lines{ii}='    var timerhandle, timeron=0;\n';
ii=ii+1; lines{ii}=['    var numberobj=' num2str(number) ';\n'];

ii=ii+1; lines{ii}='    var myLabel=new Array(numberobj);\n';

for jj=1:number
    ii=ii+1; lines{ii}=['    myLabel[' num2str(jj-1) ']="' Tag0{jj} '";\n'];
    %ii=ii+1; lines{ii}=['    myLabel[' num2str(jj-1) ']="' 'AAA' '";\n'];
end

ii=ii+1; lines{ii}='\n';
ii=ii+1; lines{ii}='    function clickalert(id)\n';
ii=ii+1; lines{ii}='    {\n';
ii=ii+1; lines{ii}='        hideall();\n';
ii=ii+1; lines{ii}='        document.getElementById("mshape" + id).setAttribute("transparency", "0.0");\n';
ii=ii+1; lines{ii}='        document.getElementById("comments").value=myLabel[id];\n';
ii=ii+1; lines{ii}='        if(timeron>0) { clearTimeout(timerhandle);} else { timeron=1; }\n';
ii=ii+1; lines{ii}='        timerhandle=setTimeout("timeron=0; showall();",1000);\n';
ii=ii+1; lines{ii}='    }\n';
ii=ii+1; lines{ii}='\n';
ii=ii+1; lines{ii}='    function hideall()\n';
ii=ii+1; lines{ii}='    {\n';
ii=ii+1; lines{ii}='        for (i=0;i<numberobj;i++)\n';
ii=ii+1; lines{ii}='        {\n';
ii=ii+1; lines{ii}='            document.getElementById("mshape" + i).setAttribute("transparency", "0.9");\n';
ii=ii+1; lines{ii}='        }\n';
ii=ii+1; lines{ii}='    }\n';
ii=ii+1; lines{ii}='    function showall()\n';
ii=ii+1; lines{ii}='    {\n';
ii=ii+1; lines{ii}='        for (i=0;i<numberobj;i++)\n';
ii=ii+1; lines{ii}='        {\n';
ii=ii+1; lines{ii}='            document.getElementById("mshape" + i).setAttribute("transparency", "0.0");\n';
ii=ii+1; lines{ii}='        }\n';
ii=ii+1; lines{ii}='    }\n';
ii=ii+1; lines{ii}=']]>\n';
str=sprintf([lines{:}]);

end
