function genmagstr(obj, varargin)
% generates magnetic structure
%
% GENMAGSTR(obj, 'option1', value1 ...)
%
% There are several ways to generate magnetic structure. The selected
% method depends on the 'mode' option, see below. The magnetic structure is
% stored in the obj.mag_str field.
%
% Input:
%
% obj       spinw class object.
%
% Options:
%
% mode       Mode how the magnetic structure is generated.
%
%            'random'
%                   generates random spin directions with zero k.
%
%            'direct'
%                   direct input of the magnetic structure from the S, k, n
%                   parameters.
%
%            'extend'
%                   Simply extend the existing or input structure (param.S),
%                   if no structure exists, random structure is generated.
%                   If defined param.S is used as starting structure for
%                   extension. If the starting structure is already
%                   extended with other size, the spins in the original
%                   unit cell are used. Magnetic ordering wavevector k will
%                   be set to zero. To generate structure with non-zero k,
%                   use 'helical' or 'direct' option. (default)
%
%            'helical'
%                   generates helical structure, starting structure
%                   is defined by param.S, the normal vector of rotation is
%                   stored in param.n, the ordering wavevector is stored in
%                   param.k. If param.S is complex, it is used as basis
%                   vectors to generate magnetic structure according to the
%                   following formula (param.n value then omitted):
%
%                   M_i(r) = R(2*pi*km*r)*M_i.
%
%                   where M_i vectors are defined in param.S. It has to
%                   contain either 1 spin direction or as many as the
%                   number of magnetic atoms in the crystallographic unit
%                   cell. In the first case, the r position is the atomic
%                   position, in the second case r is the lattice
%                   translation vector of the crystallographic cell where
%                   the moment directions are calculated.
%
%            'rotate'
%                   uniform rotation of all magnetic moments with a
%                   param.phi angle around the param.n vector. If
%                   param.phi=0, all moments are rotated so, that the first
%                   moment is parallel to param.n vector in case of
%                   collinear structure or in case of planar structure
%                   param.n defines the normal of the plane of the magnetic
%                   moments.
%
%            'func'
%                   function that defines the magnetic moments, ordering
%                   wave vector and normal vector from parameters in the
%                   following form:
%                       [M, k, n] = @(x)func(M0, x)
%                   where M is matrix with dimensions of [3 nMagExt]. k is
%                   the ordering wave vector, its size is (1,3). n is the
%                   normal vector of the spin rotation plane, its
%                   dimensions are [1 3]. The default is @gm_spherical3d.
%                   For planar magnetic structure use @gm_planar. Only
%                   param.func and param.x have to be defined for this
%                   mode.
%            'fourier'
%                   generate magnetic structure from Fourier components,
%                   usefull for multi-k structures. It creates a magnetic
%                   supercell that incorporates the magnetic structure
%                   using the following equation:
%
%                   S_{l,j} = sum_k F(k)_j*exp(-i*k*l)
%
%                   Where the summation runs over all given k-vectors, 
%                   S_{l,j} is the direction of the ordered moment of the
%                   l-th unit cell, j-th atom. The size of the generated
%                   supercell is determined by the 'nExt' option. The 'Fk'
%                   option gives the Furier components and the k-vectors in
%                   a cell in the following structure:
%                   {Fk1 k1 Fk2 k2 ...}
%                   The Fk1, Fk2 etc are complex matrices that contain the
%                   Fourier compoents on every magnetic atom in the
%                   crystallographic unit cell. They have a dimension of 
%                   [3 nMagAtom]. The k1, k2 etc are k-vectors of the
%                   Fourier compoents, with dimensions of [1 3]. Since the
%                   generated magnetic structures have to be real, the -k
%                   compoents are automatically added: F(-k) = conj(F(k)).
%                   Example input: 
%                   {[1 -1;i1 -i1;0 0] [1/2 0 0] [1 -1;0 0; i1 -i1] [0 1/2 0]}
%                   This gives a double k structure for a lattice with two
%                   magnetic atoms in the unit cell.The Fourier compoents
%                   are by default in the xyz coordinate system but if
%                   'unitS' is set to 'lu', than the moment components are
%                   assumed to be in lattice units.
% phi       Angle of rotation of the magnetic moments in rad. Default
%           value is 0.
% nExt      Number of unit cell to extend the magnetic structure,
%           dimensions are [1 3]. Default value is stored in obj.
%           If nExt is a single number, then the size of the extended unit
%           cell is automatically determined from the magnetic ordering
%           wavevector. If nExt = 0.01, then the number of unit cells is
%           determined so, that in the extended unit cell, the magnetic
%           ordering wave vector is [0 0 0], within the given 0.01 error.
% k         Magnetic ordering wavevector in r.l.u., dimensions are [1 3].
%           Default value is defined in obj.
% n         Normal vector to the spin rotation plane, dimensions are [1 3].
%           Default value [0 0 1].
% S         Direct input of the spin values, dimensions are [3 nSpin].
%           Every column defines the S_x, S_y and S_z components of the
%           moment in the xyz Descartes coodinate system.
%           Default value is stored in obj.
% Fk        Fourier compoents for a multi-k magnetic structure in a cell.
%           For description, see the 'fourier' mode description above.
%           No default value, it has to be defined if 'mode' is 'fourier'.
% unitS     Units for S, default is 'xyz', optionally 'lu' can be used,
%           in this case the input spin components are assumed to be in
%           lattice units and they will be converted to the xyz coordinate
%           system.
% epsilon   The smalles value of incommensurability that is
%           tolerated without warning. Default is 1e-5.
% func      Function that produce the magnetic moments, ordering wave
%           vector and normal vector from the param.x
%           parameters in the following form:
%             [M, k, n] = @(x)func(M0,x)
%           where M is (3,nMagExt) size matrix. k is the ordering
%           wave vector, its size is (1,3). n is the normal vector
%           of the spin rotation plane, its dimensions are [1 3]. The
%           default is @gm_spherical3d. For planar magnetic structure
%           use @gm_planar.
% x0        Input parameters for param.func function, dimensions are
%           [1 nx].
% norm      Set the length of the generated magnetic moments to be equal to
%           the spin of the magnetic atoms. Default is true.
% r0        If true and only a single spin direction is given, the spin
%           phases are determined by atomic position times k-vector, while
%           if it is false, the first spin will have zero phase.
%
% Output:
%
% The obj.mag_str field will contain the new magnetic structure.
%
% Example:
%
% USb = spinw;
% USb.genlattice('lat_const',[6.203 6.203 6.203],'angled',[90 90 90],'sym','F m -3 m')
% USb.addatom('r',[0 0 0],'S',1)
% FQ = {[0;0;1+1i] [0 0 1] [0;1+1i;0] [0 1 0] [1+1i;0;0] [1 0 0]};
% USb.genmagstr('mode','fourier','Fk',FQ,'nExt',[1 1 1])
% plot(USb,'range',[1 1 1])
%
% The above example creates the multi-q magnetic structure of USb with the
% FQ Fourier components and plots the magnetic structure.
%
% See also SPINW, SPINW.ANNEAL, SPINW.OPTMAGSTR, GM_SPHERICAL3D, GM_PLANAR.
%

if ~any(obj.atom.mag)
    error('sw:genmagstr:NoMagAtom','There is no magnetic atom in the unit cell with S>0!');
end

inpForm.fname  = {'mode'   'nExt'            'k'           'n'      };
inpForm.defval = {'extend' obj.mag_str.nExt obj.mag_str.k  [0 0 1]  };
inpForm.size   = {[1 -1]   [1 -4]            [1 3]         [1 3]    };
inpForm.soft   = {false    false             false         false    };

inpForm.fname  = [inpForm.fname  {'func'          'x0'   'norm' 'r0' }];
inpForm.defval = [inpForm.defval {@gm_spherical3d []     true   true }];
inpForm.size   = [inpForm.size   {[1 1]           [1 -3] [1 1]  [1 1]}];
inpForm.soft   = [inpForm.soft   {false           true   false  false}];

inpForm.fname  = [inpForm.fname  {'S'     'phi' 'epsilon' 'unitS' 'Fk'      }];
inpForm.defval = [inpForm.defval {[]      0     1e-5      'xyz'   cell(1,0) }];
inpForm.size   = [inpForm.size   {[-6 -2] [1 1] [1 1]     [1 -5]  [1 -7]    }];
inpForm.soft   = [inpForm.soft   {true    false false     false   false     }];

param = sw_readparam(inpForm, varargin{:});

% make string lower case
param.mode = lower(param.mode);

if isempty(param.S)
    param.S = obj.mag_str.S;
else
    % number of rows in .S has to be three
    if numel(param.S) == 3
        param.S = param.S(:);
    elseif (size(param.S,1) ~= 3)
        error('sw:genmagstr:WrongInput','Parameter S has to have dimensions of [3 nSpin]!');
    end
    
    switch lower(param.unitS)
        case 'lu'
            % convert the moments from lattice units to xyz
            param.S = obj.basisvector(true)*param.S;
            
            % convert the Fourier components if they are given
            if ~isempty(param.Fk)
                for ii = 1:2:numel(param.Fk)
                    param.Fk{ii} = obj.basisvector(true)*param.Fk{ii};
                end
            end
        case 'xyz'
            % do nothing
        otherwise
            error('sw:genmagstr:WrongInput','Parameter unitS has to be either ''xyz'' or ''lu''!');
    end
end

nExt     = double(param.nExt);

% automatic determination of the size of the extended unit cell
% if nExt is a single number
if numel(nExt) == 1
    [~, nExt] = rat(param.k,nExt);
end

mAtom    = obj.matom;
nMagAtom = size(mAtom.r,2);
nMagExt  = nMagAtom*prod(nExt);

if nMagAtom==0
    error('sw:genmagstr:NoMagAtom','There is no magnetic atom in the unit cell with S>0!');
end

% Create mAtom.Sext matrix.
mAtom    = sw_extendlattice(nExt, mAtom);

% Magnetic ordering wavevector
k = param.k;

% Axis of rotation, size (1,3)
n = (param.n)/norm(param.n);

% If the magnetic structure is not initialized start with a random one.
if strcmp(param.mode,'extend') && (nMagAtom > size(param.S,2))
    param.mode = 'random';
    % this warning is not necessary
    %warning('sw:genmagstr:WrongInitialStructure','No magnetic structure is defined, random structure is created instead!')
end

if obj.symb
    k = sym(k);
    param.S = sym(param.S);
    n = sym(n);
end

switch param.mode
    case 'extend'
        % Extend the unit cell if:
        % -the new number of extended cells does not equal to the number of
        %  cells defined in obj
        % -the number of spins stored in obj is not equal to the number
        %  of spins in the final structure
        if any(obj.mag_str.N_ext - int32(param.nExt)) || (size(param.S,2) ~= nMagExt)
            S = param.S(:,1:nMagAtom);
            S = repmat(S,[1 prod(nExt)]);
        else
            S = param.S;
        end
        k = [0 0 0];
    case 'random'
        % Create random spin directions.
        S = randn(nMagExt,3);
        S = bsxfun(@rdivide,S,sqrt(sum(S.^2,2)));
        S = bsxfunsym(@times,S,mAtom.Sext')';
        k = [0 0 0];
    case 'helical'
        S0 = param.S;
        % Magnetic ordering wavevector in the extended unit cell.
        kExt = k.*nExt;
        % Warns about the non sufficient extension of the unit cell.
        % we substitute random values for symbolic km
        skExt = sw_sub1(kExt,'rand');
        if any(abs(skExt-round(skExt))>param.epsilon) && prod(nExt) > 1
            warning('sw:genmagstr:UCExtNonSuff','In the extended unit cell k is still larger than epsilon!');
        end
        % Number of spins in the input.
        nSpin = size(param.S,2);
        
        if (nSpin~= nMagAtom) && (nSpin==1)
            % Single defined spin, use the atomic position.
            if param.r0
                r = mAtom.RRext;
            else
                r = bsxfun(@minus,mAtom.RRext,mAtom.RRext(:,1));
            end
        elseif nSpin == nMagAtom
            % First crystallographic unit cell defined, use only unit cell
            % position.
            r = bsxfun(@rdivide,floor(bsxfun(@times,mAtom.RRext,nExt')),nExt');
        else
            error('sw:genmagstr:WrongNumberSpin','Wrong number of input spins!');
        end
        % Angles of rotation for each unit cell.
        phi = kExt*r*2*pi;
        
        % Spin in the extended unit cell.
        S = zeros(3,nMagExt);
        if obj.symbolic
            S = sym(S);
        end
        
        if isreal(param.S) || isa(param.S,'sym')
            % Rotate spins for each unit cell.
            for ii = 1:nMagExt
                selS    = S0(:,mod(ii-1,nSpin)+1);
                S(:,ii) = sw_rot(n,phi(ii),selS);
            end
        else
            error('sw:genmagstr:WrongInput','For complex Fourier components use the ''mode'' ''fourier'' option!');
        end
    case 'direct'
        % Direct input of the magnetic moments.
        S = param.S;
        if size(S,2)==nMagAtom
            % single unit cell
            nExt = [1 1 1];
        end
        
        if (size(S,1) ~= 3) || (size(S,2) ~= nMagExt)
            error('sw:genmagstr:WrongSpinSize','Wrong size of param.S!');
        end
    case 'rotate'
        S   = param.S;
        if param.phi == 0
            % The starting vector, size (1,3):
            S1 = sw_nvect(S);
            % Axis of rotation defined by the spin direction
            nRot  = cross(n,S1);
            % Angle of rotation.
            phi = -atan2(norm(cross(S1,n)),dot(S1,n));
        else
            nRot = n;
            % Angle of rotation.
            phi = param.phi(1);
        end
        % Rotate the spins.
        S = sw_rot(nRot,phi,S);
        n = obj.mag_str.n;
        k = obj.mag_str.k;
    case 'func'
        S = mAtom.S;
        S = repmat(S,[prod(nExt) 1]);
        
        if obj.symbolic
            [S, k, n] = param.func(sym(S), sym(param.x0));
        else
            [S, k, n] = param.func(S,param.x0);
        end
    case 'fourier'
        % generate supercell from Fourier components
        % keeps the final k-vector zero
        Fk = param.Fk;
        if isempty(Fk) || ~iscell(Fk)
            error('sw:genmagstr:WrongInput','Wrong ''Fk'' option that defines the Fourier components!');
        end
        
        % number of moments for the Fourier components are defined
        nFourier = size(Fk{1},2);
        nQ = numel(Fk)/2;
        
        if (nFourier~= nMagAtom) && (nFourier==1)
            % Single defined moment, use the atomic position in l.u.
            RR = bsxfun(@times,mAtom.RRext,nExt');
        elseif nFourier == nMagAtom
            % First crystallographic unit cell defined, use only unit cell
            % position in l.u.
            RR = floor(bsxfun(@times,mAtom.RRext,nExt'));
        else
            error('sw:genmagstr:WrongNumberComponent','Wrong number of input Fourier components!');
        end
        
        % no moments
        S = RR*0;
        % number of cells in the supercell
        nCell = prod(nExt);
        
        % multiply the Fourier components with the spin quantum number
        % TODO
        
        
        for ii = 1:2:(2*nQ)
            % F(k)
            S = S + bsxfunsym(@times,repmat(Fk{ii},[1 nCell*nMagAtom/nFourier]),exp(1i*Fk{ii+1}*RR*2*pi))/2;
            % conj(F(k))
            S = S + bsxfunsym(@times,repmat(conj(Fk{ii}),[1 nCell*nMagAtom/nFourier]),exp(-1i*Fk{ii+1}*RR*2*pi))/2;
            
        end
        S = real(S);

        k = [0 0 0];
        
        warning('sw:genmagstr:Approximation',['The generated magnetic '...
            'structure is only an approximation of the input multi-q'...
            ' structure on a supercell!'])
        n = [0 0 1];
        
    otherwise
        error('sw:genmagstr:WrongMode','Wrong param.mode value!');
end

% normalize the magnetic moments
if param.norm
    normS = sqrt(sum(S.^2,1))./repmat(mAtom.S,[1 prod(nExt)]);
    normS(normS==0) = 1;
    S = bsxfunsym(@rdivide,S,normS);
    %if any(isnan(S(:)))
    %    error('spinw:genmagstr:WrongMoments','Zero magnetic moments cannot be normalized!');
    %end
    
end

% simplify expressions
if obj.symbolic
    S = simplify(sym(S));
    k = simplify(sym(k));
    n = simplify(sym(n));
end

mag_str.N_ext = int32(nExt(:))';
mag_str.k     = k(:)';
mag_str.S     = S;
mag_str.n     = n(:)';

obj.mag_str   = mag_str;

end