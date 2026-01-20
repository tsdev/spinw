function  Eout = energy(obj, varargin)
% calculates the ground state energy
%
% ### Syntax
%
% `E = energy(obj,Name,Value)`
%
% ### Description
%
% `E = energy(obj,Name,Value)` calculates the classical ground state energy
% per spin. The calculation correctly takes into account the magnetic
% supercell. The function gives correct results on single-k magnetic
% structures even defined on magnetic supercells. For multi-k magnetic
% structures first a definition of a larger supercell is necessary where an
% effective $k=0$ representation is possible.
%
% ### Examples
%
% After optimising the magnetic structure (by minimizing the ground state
% energy), the energy per spin is calculated. This can be compared to
% different ground state structures to decide which is the right classical
% ground state of the magnetic model in cryst. Here we use the triangular
% lattice antiferromagnet where we define the magnetic structure on a
% $[3\times 3]$ magnetic supercell where the optimal structure (120\\deg
% angle between neighboring spins) has a 0 propagation vector. In this case
% the exact energy is $3\cdot 1^2\cdot \cos(120^\circ) = -1.5$.
%
% ```
% >>cryst = sw_model('triAF',1)
% >>cryst.genmagstr('mode','random','nExt',[3 3 1])
% >>cryst.optmagsteep('nRun',10)
% >>cryst.energy>>
% ```
%
% ### Input Arguments
%
% `obj`
% : [spinw] object.
%
% ### Name-Value Pair Arguments
%
% `'epsilon'`
% : The smallest value of incommensurability that is tolerated
%   without warning. Default is $10^{-5}$.
%
% ### Output Arguments
%
% `E`
% : Energy per moment (anisotropy + exchange + Zeeman energy).
%
% {{warning The calculated energy can be wrong for incommensurate
% structures. For example a structure where the spins are rotating in $XY$
% plane with an incommensurate wavevector of $(1/3,0,0)$. The function only
% calculates the anisotropy energy in the first unit cell, that is for
% single spin $E_{aniso} = A_{xx}\cdot S_{xx}^2+A_{yy}\cdot S_{yy}^2$.
% While the anisotropy energy in reality is independent of the spin
% orientation in the $XY$ plane $E_{aniso}=3S\cdot (A_{xx}+A_{yy})/2$. Thus
% using `spinw.energy` on incommensurate structures together with single
% ion anisotropy one has to be carefull! In the triangular case one has to
% extend the unit cell to `nExt = [3 3 1]` (in the hexagonal setting), in
% this case the energy will be correct.}}
%
% ### See Also
%
% [spinw] \| [spinw.anneal] \| [spinw.newcell]
%

inpForm.fname  = {'epsilon'};
inpForm.defval = {1e-5     };
inpForm.size   = {[1 1]    };

param = sw_readparam(inpForm,varargin{:});

[SS, SI] = obj.intmatrix;

% add dipolar terms
SS.all = [SS.all SS.dip];

magStr = obj.magstr;
nExt   = magStr.N_ext;
kExt   = magStr.k.*nExt;

% Incommensurate structure in the extended unit cell.
if obj.symbolic
    kTest = sw_sub1(kExt,0.452524243);
else
    kTest = kExt;
end
incomm = any(abs(kTest-round(kTest))>param.epsilon);

khalf = ~any(mod(kTest*2,1));

M0      = magStr.S;
nMagExt = size(M0,2);

% Anisotropy energy can be wrong
isAniso = any(sw_sub1(SI.aniso(:)));
gMat    = sw_sub1(SI.g);
isGoff  = any(any(any(abs(gMat-gMat(1)*repmat(eye(3),[1 1 nMagExt]))>param.epsilon)));

if  incomm && (isAniso || any(SI.field) || isGoff) && ~isreal(obj.mag_str.F) && ~khalf
    warning('spinw:energy:AnisoFieldIncomm','The calculated energy might be wrong due to incommensurate structure!');
end

if nMagExt>0
    
    dR    = [SS.all(1:3,:) zeros(3,nMagExt)];
    atom1 = [SS.all(4,:)   1:nMagExt];
    atom2 = [SS.all(5,:)   1:nMagExt];
    JJ    = cat(3,reshape(SS.all(6:14,:),3,3,[]),SI.aniso);
    
    M1 = M0(:,atom1);
    M2 = M0(:,atom2);
    M2cmplx = obj.mag_str.F(:,atom2);
    
    % For incommensurate structures in the extended unit cell rotates M2
    % moments.
    if incomm && any(any(dR))
        %n = magStr.n;
        dRIdx = find(any(dR));
        for ii = 1:numel(dRIdx)
            %M2(:,dRIdx(ii)) = sw_rot(n, kExt*dR(:,dRIdx(ii))*2*pi, M2(:,dRIdx(ii)));
            % TODO check - sign in front of phase
            M2(:,dRIdx(ii)) = real(bsxfun(@times,M2cmplx(:,dRIdx(ii)),exp(-1i*kExt*dR(:,dRIdx(ii))*2*pi)));
        end
    end
    
    Ml = repmat(permute(M1,[1 3 2]),[1 3 1]);
    Mr = repmat(permute(M2,[3 1 2]),[3 1 1]);

    Q = Ml.*JJ.*Mr;
    if isa(Ml,'sym') || isa(Ml,'sym') || isa(Ml,'sym')
        [n1, n2, n3] = size(Q);
        sum_3 = 0;
        for k = 1:n3
            sum_3 = sum_3 + Q(:,:,k);
        end
    else
        sum_3 = sum(Q,3);
    end
    exchE   =  sum(sum(sum_3,2),1);

    % TODO
    % correct energy for twins
    ZeemanE = -sum(SI.field*permute(mmat(SI.g,permute(M0,[1 3 2])),[1 3 2]),2)*obj.unit.muB;
    
    % energy per spin
    E = (exchE + ZeemanE)/nMagExt;
    
else
    % for undefined magnetic structure, energy returns NaN
    E = NaN;
end

if nargout > 0
    Eout = E;
else
    if obj.symbolic
        fprintf('Ground state energy (1/spin):\n')
        pretty(sym('E')==E);
    else
        fprintf(['Ground state energy: %5.3f ' obj.unit.label{2} '/spin.\n'],E);
    end
end

end