function op = symop(obj)
% generates the bond symmetry operators
% 
% ### Syntax
% 
% `op = symop(obj)`
% 
% ### Description
% 
% `op = symop(obj)` generates the rotation matrices that transform single
% ion anisotropy, g-tensor and exchange interaction matrices between
% symmetry equivalent positions (on atoms or bond centers). The results are
% cached.
% 
% ### See Also
% 
% [spinw.intmatrix]
%

% bonds
coupling = obj.coupling;

% bond list
bondList = double([coupling.dl; coupling.atom1; coupling.atom2; coupling.idx]);

% find the last symmetry generated matrix
lastSym = find(coupling.idx <= coupling.nsym,1,'last');

if isempty(obj.cache.symop)
    % generate the symmetry operators and fill the cache
    if coupling.nsym > 0
        % transformation matrix between l.u. and xyz coordinate systems
        A = obj.basisvector(false,obj.symbolic);
        
        
        % generate symmetry operators for anisotropy matrice using the space group symmetry
        % generate the rotation matrices
        if obj.symbolic
            [~, aIdx, opInfo] = swsym.position(obj.lattice.sym,obj.unit_cell.r(:,~sw_always(obj.unit_cell.S==0)));
        else
            [~, aIdx, opInfo] = swsym.position(obj.lattice.sym,obj.unit_cell.r(:,obj.unit_cell.S>0));
        end
        
        rotOpA = opInfo.opmove;
        % find the first atoms on every orbit
        atomSel = [true logical(diff(aIdx))];
        % take the atoms that generate the rest of the atoms on orbit
        firstAtom = aIdx(atomSel);
        % create unit operators on the first atoms
        for ii = 1:numel(firstAtom)
            % operator on the first bond
            Op0 = rotOpA(:,:,find(aIdx==ii,1));
            % correct
            rotOpA(:,:,aIdx == ii) = mmat(rotOpA(:,:,aIdx == ii),inv(Op0));
        end
        
        % convert rotation operators to xyz Cartesian coordinate system
        rotOpA = mmat(A,mmat(rotOpA,inv(A)));
        
        % generate symmetry operators for exchange matrices
        % first positions of the couplings with identical idx values used to
        % generate the coupling matrices for the rest
        % only calculate for the symmetry generated bonds
        bondSel = [true logical(diff(bondList(6,:)))] & bondList(6,:)<= coupling.nsym;
        % keep the bonds that will generate the space group operators
        firstBond = bondList(1:5,bondSel);
        rotOpB = zeros(3,3,lastSym);
        % select rotation matrices for each generated coupling
        bIdx = 0;
        for ii = 1:size(firstBond,2)
            [~, rotIdx] = swsym.bond(obj.matom.r,obj.basisvector, firstBond(:,ii), obj.lattice.sym, 1e-5);
            % operator on the first bond
            Op0 = obj.lattice.sym(:,1:3,find(rotIdx,1));
            rotOpB(:,:,bIdx+(1:sum(rotIdx))) = mmat(obj.lattice.sym(:,1:3,rotIdx),inv(Op0));
            bIdx = bIdx + sum(rotIdx);
        end
        
        % convert rotation operators to xyz Cartesian coordinate system
        rotOpB = mmat(A,mmat(rotOpB,inv(A)));
        
        % save to the cache
        obj.cache.symop.sion = rotOpA;
        obj.cache.symop.bond = rotOpB;
    else
        % no symmetry
        obj.cache.symop.sion = zeros(3,3,0);
        obj.cache.symop.bond = zeros(3,3,0);
    end
    % add listener to lattice and unit_cell fields
    obj.addlistenermulti(2);
end

op = obj.cache.symop;

end