function [xLabel, xAxis] = sw_label(spectra)
% [xLabel, xAxis] = SW_LABEL(spectra) returns the label and axis vector for
% the x-axis for momentum transfer scans linear in reciprocal space.
%

hkl = spectra.hkl';
hklA = spectra.hklA';

nQ  = size(hkl,1);

% distance between start and end points
dk0 = hkl(1,:) - hkl(end,:);

% determine whether it is line scan
if abs(dk0(1))>1e-5
    dk    = dk0/dk0(1);
    xAxis = hkl(:,1);
elseif abs(dk0(2))>1e-5
    dk    = dk0/dk0(2);
    xAxis = hkl(:,2);
elseif abs(dk0(3))>1e-5
    dk    = dk0/dk0(3);
    xAxis = hkl(:,3);
end

% parse curve into straight lines
qStep = hkl(2:end,:)-hkl(1:(end-1),:);
qStep = bsxfun(@rdivide,qStep,sqrt(sum(qStep.^2,2)));
qCurve = sum(qStep(2:end,:).*qStep(1:end-1,:),2);
qIdx = find(qCurve<0.97)+1;

if numel(qIdx) == 0
    linescan = 0;
elseif numel(qIdx)/nQ<0.2
    % line scan with straight pieces
    linescan = 1;
else
    % curved scan
    linescan = 2;
end

switch linescan
    case 0
        % single linear scan in r.l.u.
        k0 = hkl(1,:);
        
        changeX = false;
        inA = sqrt(sum(((spectra.hklA(:,1) - spectra.hklA(:,end)).^2)))/abs(xAxis(end)-xAxis(1));
        xiLabel = cell(1,3);
        for ii = 1:3
            if abs(k0(ii)) > 1e-3
                if abs(dk(ii)) > 1e-3
                    if abs(dk(ii)-1)<1e-3
                        xiLabel{ii} = sprintf('%.4g+\\xi',k0(ii));
                        if ~changeX
                            xAxis = xAxis - k0(ii);
                            changeX = true;
                        end
                    elseif abs(dk(ii)+1)<1e-3
                        xiLabel{ii} = sprintf('%.4g-\\xi',k0(ii));
                        if ~changeX
                            xAxis = xAxis + k0(ii);
                            changeX = true;
                        end
                        
                    else
                        xiLabel{ii} = sprintf('%.4g%+.4g\\xi',k0(ii),dk(ii));
                        if ~changeX
                            xAxis = xAxis - k0(ii);
                            changeX = true;
                        end
                        
                    end
                    
                else
                    xiLabel{ii} = sprintf('%.4g',k0(ii));
                end
            else
                if abs(dk(ii)) > 1e-3
                    if abs(dk(ii)-1)<1e-3
                        xiLabel{ii} = '\xi';
                    elseif abs(dk(ii)+1)<1e-3
                        xiLabel{ii} = '-\xi';
                    else
                        xiLabel{ii} = sprintf('%.4g\\xi',dk(ii));
                    end
                else
                    xiLabel{ii} = '0';
                end
            end
        end
        if size(hkl,1) == 1
            xLabel = ['[' xiLabel{1} ',' xiLabel{2} ',' xiLabel{3} ']'];
        else
            xLabel = ['[' xiLabel{1} ',' xiLabel{2} ',' xiLabel{3} '] in ' sprintf('%.5g',inA) ' A^{-1}'];
        end
        
    case 1
        % use inverse Angstrom for the x-axis scaling
        qIdx = [1;qIdx;nQ];
        xAxis = 0;
        for ii = 2:length(qIdx)
            hklAdist = sqrt(sum((hklA(qIdx(ii),:)-hklA(qIdx(ii-1),:)).^2));
            qAdd = linspace(0,hklAdist,qIdx(ii)-qIdx(ii-1)+1);
            qAdd = qAdd(2:end);
            xAxis = [xAxis qAdd+xAxis(end)]; %#ok<AGROW>
        end
        % create labels for line pieces
        xLabel = cell(1,length(qIdx));
        for ii = 1:length(qIdx)
            hkl(abs(hkl)<1e-3) = 0;
            xLabel{ii} = ['(' num2str(hkl(qIdx(ii),1)) ',' num2str(hkl(qIdx(ii),2)) ',' num2str(hkl(qIdx(ii),3)) ')'];
        end
        xLabel{end+1} = xAxis(qIdx);
        1;
    case 2
        xLabel = 'Momentum transfer';
        xAxis  = linspace(0,1,nQ);
end