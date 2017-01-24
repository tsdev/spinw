function sw_status(percent,varargin)
% timer function that displays also the remaining time
%
% SW_STATUS(percent, {mode},{fid},{title})
%
% Input:
%
% percent   Percentage of the calculation that is done.
% mode      Determines the time estimation, optional parameter:
%               1   Starts the time estimation.
%               0   Displays of the remaining time. (default)
%               2   Calculation finished.
% fid       File identifier to print the output:
%               1   Text output to the Command Window. Default.
%               2   Graphical output, using the waitbar() function.
%
% See also WAITBAR.
%

if nargin == 0
    help sw_status
    return
end

if nargin > 2 && ~isempty(varargin{2})
    fid = varargin{2};
else
    fid = swpref.getpref('tid',[]);
end

if nargin>3
    title0 = varargin{3};
else
    title0 = 'sw_status';
end

if ~ismember(fid,[1 2])
    return
end

if nargin > 1
    start = varargin{1};
else
    start = 0;
end

switch start
    case 1
        % start the time estimation
        tic
        switch fid
            case 1
                fprintf([repmat(' ',[1 40]) '\n']);
            case 2
                hBar = waitbar(0,'sw_status() initializing...');
                hBar.HandleVisibility='on';
                hBar.Tag = 'sw_status';
                hBar.Name = title0;
        end
    case 0
        % refresh the displayed time
        etime = double(toc);
        if percent == 0
            percent = 1e-5;
        end
        rtime = (100-percent)/percent*etime;
        hou = floor(rtime/60^2);
        rtime = rtime-hou*60^2;
        min = floor(rtime/60);
        sec = floor(rtime - min*60);
        switch fid
            case 1
                fprintf([repmat('\b',[1 41]) '%6.2f%%, remained: %03d:%02d:%02d (HH:MM:SS).\n'],...
                    percent,hou,min,sec);
            case 2
                hBar = findobj('Tag','sw_status');
                if ~isempty(hBar)
                    waitbar(percent/100,hBar,sprintf('%6.2f%%, remained: %03d:%02d:%02d (HH:MM:SS)',percent,hou,min,sec))
                end
        end
        
    case  2
        % finish extimation
        etime = double(toc);
        hou = floor(etime/60^2);
        etime = etime-hou*60^2;
        min = floor(etime/60);
        etime = etime-min*60;
        sec = floor(etime);
        %tho = floor((etime-sec)*1000);
        %fprintf('Finished in %02d:%02d:%02d.%03d (HH:MM:SS.FFF).\n',hou,min,sec,tho);
        switch fid
            case 1
                fprintf(repmat('\b',1,40+1));
                fprintf('Calculation is finished in %02d:%02d:%02d (hh:mm:ss).\n',hou,min,sec);
            case 2
                hBar = findobj('Tag','sw_status');
                delete(hBar);
                fprintf('Calculation is finished in %02d:%02d:%02d (hh:mm:ss).\n',hou,min,sec);
        end
end

end