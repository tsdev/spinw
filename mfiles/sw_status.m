function sw_status(percent,varargin)
% timer function that displays also the remaining time
%
% SW_STATUS(percent, {mode})
%
% Input:
%
% percent   Percentage of the calculation that is done.
% mode      Determines the time estimation, optional parameter:
%               1   Starts the time estimation.
%               0   Displays of the remaining time. (default)
%               2   Calculation finished.
%

if nargin == 0
    help sw_status;
    return;
end

if nargin > 1
    start = varargin{1};
else
    start = 0;
end

if start == 0
    etime = double(toc);
    if percent == 0
        percent = 1e-5;
    end
    rtime = (100-percent)/percent*etime;
    hou = floor(rtime/60^2);
    rtime = rtime-hou*60^2;
    min = floor(rtime/60);
    sec = floor(rtime - min*60);
    fprintf([repmat('\b',[1 40]) '%6.2f%%, remained: %02d:%02d:%02d (HH:MM:SS).\n'],...
        percent,hou,min,sec);
elseif start == 1
    tic
    fprintf([repmat(' ',[1 39]) '\n']);
elseif start == 2
    etime = double(toc);
    fprintf(repmat('\b',1,40));
    hou = floor(etime/60^2); 
    etime = etime-hou*60^2;
    min = floor(etime/60); 
    etime = etime-min*60;
    sec = floor(etime);
    %tho = floor((etime-sec)*1000);
    %fprintf('Finished in %02d:%02d:%02d.%03d (HH:MM:SS.FFF).\n',hou,min,sec,tho);
    fprintf('Finished in %02d:%02d:%02d (hh:mm:ss).\n',hou,min,sec);
end
end