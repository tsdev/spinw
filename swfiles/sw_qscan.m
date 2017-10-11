function qOut = sw_qscan(qLim)
% creates continuous line between coordinates
% 
% ### Syntax
% 
% `qOut = sw_qscan(qLim)`
% 
% ### Description
% 
%  `qOut = sw_qscan(qLim)` generates connected lines between given
%  positions in $n$-dimensional space ($n>1$). The output contains equally
%  spaced points along each line section in a matrix, by default 100
%  points. The function can be used to generates points along a path
%  defined by corner points.
%
% ### Input Arguments
%
% `qLim`
% : Cell that contains row vectors with $n$ elements each and optionally an
%   additional integer, e.g. `{[0 0] [1 0] 201}`.
%
% ### Examples
% 
% To generate a path in the Brillouin-zone between the $(0,0,0)$, $(1,0,0)$
% and $(1,1,0)$ points with 501 points per line segment use:
%
% ```
% >>Q = sw_qscan({[0 0 0] [1 0 0] [1 1 0] [0 0 0] 501});
% >>>figure
% >>plot(Q(1,:),Q(2,:),'linewidth',2)
% >>xlabel H
% >>ylabel K
% >>>axis([-1 2 -1 2])
% >>>grid on
% >>snapnow
% ```
%
% ### See Also
%
% [sw_qgrid]
%

if nargin == 0
    help sw_qscan
    return
end

if ~iscell(qLim)
    if size(qLim,1) > 3 || size(qLim,3)>1
        error('sw_qscan:WrongInput','The dimensions of the q-vector list are wrong!')
    end
    qOut = qLim;
    return
end

if numel(qLim{end}) == 1
    nQ = qLim{end};
    qLim = qLim(1:end-1);
else
    nQ = 100;
end

if iscell(qLim) && length(qLim)>1
    qOut = zeros(length(qLim{1}),0);
    for ii = 2:length(qLim)
        q1 = reshape(qLim{ii-1},[],1);
        q2 = reshape(qLim{ii},  [],1);
        
        if nQ > 1
            qOut = [qOut bsxfun(@plus,bsxfun(@times,q2-q1,linspace(0,1,nQ)),q1)]; %#ok<AGROW>
        else
            qOut = (q2+q1)/2;
        end
        if ii<length(qLim)
            qOut = qOut(:,1:end-1);
        end
    end
elseif iscell(qLim)
    qOut = qLim{1};
else
    qOut = qLim;
end

end