function sw_deploytest(command)
% code to compile that runs with pymatbridge
%
% to compile:
% mcc -m -d deploytest -o sw_deploytest sw_deploytest.m

eval(command)

end