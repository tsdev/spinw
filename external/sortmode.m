function [Dseq,Vseq] = sortmode(Dseq,Vseq,deg)
% sort eigenvalues and eigenvectors
%
% ### Syntax
%
% `[Dseq,Vseq] = sortmode(Dseq,Vseq)`
%
% `[Dseq,Vseq] = sortmode(___,deg)`
%
% ### Description
%
% `[Dseq,Vseq] = sortmode(Dseq,Vseq)` sorts the eigenvalues and
% eigenvectors stored in `Dseq` and `Vseq` matrices. The code tries to find
% continuous dispersion lines first based on smoothness (by fitting
% dispersion up to 4th degree polynomial) and at the crossing points of
% multiple dispersion lines it tries to minimize the variation of the
% eigenvectors.
%
% {{note If the dispersion is extensively degenerate, the modes are sorted
% only according to the smoothness of the dispersion.}}
%
% {{warning The code assumes that the input dispersion matrix is sorted in
% energy.}}
%
% `[Dseq,Vseq] = sortmode(___,deg)` control wether the spectrum is
% extensively degenerate or not by setting `deg` `true` or `false`
% respectively.
%
% ### Input Arguments
%
% `Desq`
% : Dispersion stored in a matrix with dimensions of $[n_{mode}\times
%   n_{point}]$.
%
% `Vseq`
% : Eigenvectors or any other value corresponding to the eigenvalues stored
%   in matrix with dimensions of $[n_{vec}\times n_{mode}\times
%   n_{point}]$, where $v_{vec}=n_{mode}$ if `Vseq` stores the
%   eigenvectors, otherwise it can be of arbitrary size.
%
% ### Output Arguments
%
% `Dseq`, `Vseq`
% : Same as the input matrices except that they are sorted along the
%   $n_{mode}$ dimension.
%

% number of matrices to sort
N          = size(Dseq,2);

if N<3
    % nothing to sort
    return
end

nMode      = size(Dseq,1);
idx        = zeros(nMode,N);
idx(:,1:2) = repmat((1:nMode)',[1 2]);

if nargin<3
    % determine the level of degeneracy in the spectrum
    deg = sum(abs(diff(sort(real(Dseq(:)))))./max(real(Dseq(:)))<1e-5)/numel(Dseq);
    if deg > 5e-2
        deg = true;
    else
        deg = false;
    end
end

% calculate polyN extrapolation
%poldeg = @(N)subs(poly2sym(linsolve(sym(fliplr(((0:N)').^(0:N))),sym('y',[1 N+1],'real')')),sym('x'),N+1);

for ii = 3:N
    % new energy values
    if ii == 3
        Ep  = 2*Dseq(:,ii-1)-Dseq(:,ii-2);
    else
        % elseif ii == 4
        Ep  = 3*Dseq(:,ii-1)-3*Dseq(:,ii-2)+Dseq(:,ii-3);
        % elseif ii == 5
        %   Ep  = 4*Dseq(:,ii-1)-6*Dseq(:,ii-2)+4*Dseq(:,ii-3)-Dseq(:,ii-4);
        % elseif ii == 6
        %   Ep  = 5*Dseq(:,ii-1)-10*Dseq(:,ii-2)+10*Dseq(:,ii-3)-5*Dseq(:,ii-4)+Dseq(:,ii-5);
        %   Ep  = 6*Dseq(:,ii-1)-15*Dseq(:,ii-2)+20*Dseq(:,ii-3)-15*Dseq(:,ii-4)+6*Dseq(:,ii-5)-Dseq(:,ii-6);
    end
    % calculate distance matrix
    idx(:,ii) = munkres(abs(bsxfun(@minus,Ep,Dseq(:,ii)')));
    % resort modes
    Dseq(:,ii)   = Dseq(idx(:,ii),ii);
    Vseq(:,:,ii) = Vseq(:,idx(:,ii),ii);
    
    if ~deg
        % invert idx
        [~,idx(:,ii)] = sort(idx(:,ii));
        
        % find mode crossings only for non-degenerate spectrum
        modeIdx = find(idx(:,ii)-idx(:,ii-1));
        
        if ~isempty(modeIdx)
            iSel = idx(modeIdx,ii+[-1 0]);
            
            % mode indices that participate in the crossing
            [uMode,~,uiSel]  = unique(iSel(:));
            uiSel = reshape(uiSel,2,[]);
            color  = zeros(1,numel(uMode));
            
            for jj = 1:size(iSel,2)
                % go over all links
                if color(uiSel(1,jj))>0
                    color(uiSel(2,jj)) = color(uiSel(1,jj));
                elseif color(uiSel(2,jj))>0
                    color(uiSel(1,jj)) = color(uiSel(2,jj));
                else
                    color(uiSel(:,jj)) = max(color)+1;
                end
            end
            
            % apply munkres on the crossing points
            for jj = 1:max(color)
                % mode index for a selected crossing
                index0 = uMode(color==jj);
                % find eigenvectors for these modes
                selP   = permute(Vseq(:,index0,ii+[-1 0]),[2 3 1]);
                D1     = sum(abs(bsxfun(@minus,selP(:,1,:),permute(selP(:,2,:),[2 1 3]))).^2,3);
                D2     = sum(abs(bsxfun(@plus,selP(:,1,:),permute(selP(:,2,:),[2 1 3]))).^2,3);
                indexm = munkres(min(cat(3,D1,D2),[],3));
                % new modes: index0 --> index0(indexm)
                index1 = index0(indexm);
                % flip the modes
                Vseq(:,index0,ii) = Vseq(:,index1,ii);
                Dseq(index0,ii) = Dseq(index1,ii);
            end
            
        end
    end
    
end

end