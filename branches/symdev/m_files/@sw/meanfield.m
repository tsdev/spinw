function MF = meanfield(obj, varargin)
% returns different mean field values

inpForm.fname  = {'epsilon'};
inpForm.defval = {1e-5     };
inpForm.size   = {[1 1]    };

param = sw_readparam(inpForm,varargin{:});

[SS, SI] = obj.intmatrix('conjugate',true);

nExt = double(obj.mag_str.N_ext);
kExt = obj.mag_str.k.*nExt;

% Incommenasurate structure in the extended unit cell.
incomm = any(abs(kExt-round(kExt))>param.epsilon);

M0      = obj.mag_str.S;
nMagExt = size(M0,2);

% Anisotropy energy can be wrong.
if  incomm && (any(any(any((SI.aniso)))) || (any(SI.field)) || (any(any(any(abs(SI.g-SI.g(1)*repmat(eye(3),[1 1 nMagExt]))>param.epsilon)))))
    warning('sw:energy:AnisoFieldIncomm','The calculated energy might be wrong due to incommensurate structure!');
end

if nMagExt>0

    dR    = SS.all(1:3,:);
    atom1 = SS.all(4,:);
    atom2 = SS.all(5,:);
    JJ    = reshape(SS.all(6:end,:),3,3,[]);
    
    %M1 = M0(:,atom1);
    M2 = M0(:,atom2);
    
    % For incommensurate structures in the extended unit cell rotates M2
    % moments.
    if incomm && any(any(dR))
        n = obj.mag_str.n;
        dRIdx = find(any(dR));
        for ii = 1:length(dRIdx)
            M2(:,dRIdx(ii)) = sw_rot(n, kExt*dR(:,dRIdx(ii))*2*pi, M2(:,dRIdx(ii)));
        end
    end
    
    %Ml = repmat(permute(M1,[1 3 2]),[1 3 1]);
    Mr = repmat(permute(M2,[3 1 2]),[3 1 1]);
    
    exchE   =  squeeze(sum(JJ.*Mr,2));
    
    field = zeros(3,nMagExt);
    for ii = 1:3
        field(ii,:) = accumarray(atom1',exchE(ii,:),[nMagExt 1]);
    end
    
    MF.field = field;
    
else
    % for undefined magnetic structure, energy returns NaN
    MF.field = NaN;
end


end