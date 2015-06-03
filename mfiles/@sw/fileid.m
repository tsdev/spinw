function fidOut = fileid(obj,fid)
% determines where the text out is written
%
% FILEID(obj,fid)
%
% Determines the text output of all sw class methods. Default
% is 1, where all output is printed onto the MATLAB Command
% Window.
%
% fidOut = FILEID(obj)
%
% Outputs the stored fileID value.
%
% See also SW, SW.DISPLAY.
%

if nargin > 1
    obj.fid = fid;
    return
end

fidOut = obj.fid;

end