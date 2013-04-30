function qOut = sw_qscan(qLim)
% qOut = SW_QSCAN(qLim) creates linear scans between Q points in qLim.
%
% Example:
% qLim = {[0 1 0] [0 0 0]}
% If the last element of qLim is a scalar, it defines the number of point
% in each linear scan, by defult this value is 100.
% qLim = {[0 1 0] [0 0 0] 50}
%

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
        
        qOut = [qOut bsxfun(@plus,bsxfun(@times,q2-q1,linspace(0,1,nQ)),q1)]; %#ok<AGROW>
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