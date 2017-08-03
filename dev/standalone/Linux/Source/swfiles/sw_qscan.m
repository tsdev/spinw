function qOut = sw_qscan(qLim)
% creates linear scans between Q points in 3D
%
% qOut = SW_QSCAN(qLim)
%
% Example:
%
% qLim = {[0 1 0] [0 0 0]}
% If the last element of qLim is a scalar, it defines the number of point
% in each linear scan, by default this value is 100.
% qLim = {[0 1 0] [0 0 0] 50}
%

if nargin == 0
    help sw_qscan;
    return;
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