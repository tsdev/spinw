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
    rtime = (100-percent)/percent*etime;
    fprintf('\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b%6.2f%%, remained: %8.1f minute.\n',...
        mod(percent,100),rtime/60);
elseif start == 1
    tic
    fprintf('                                   \n');
elseif start == 2
    etime = double(toc);
    fprintf('\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b');
    fprintf('Finished in %8.1f min.\n',etime/60);
end
end