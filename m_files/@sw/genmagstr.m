function obj = genmagstr(obj, varargin)
% generates magnetic structure
%
% obj = GENMAGSTR(obj, 'option1', value1 ...)
%
% There are several ways to generate magnetic structure. The method depends
% on the 'mode' option, see below. The magnetic structure is stored in the
% sw object mag_str property.
%
% Options:
%
% mode       Mode how the magnetic structure is generated.
%
%                random:
%                   generates random spin directions with zero k.
%
%                direct:
%                   direct input of the magnetic structure from the S, k, n
%                   parameters.
%
%                extend:
%                   Simply extend the existing or input structure (param.S),
%                   if no structure exists, random structure is generated.
%                   If defined param.S is used as starting structure for
%                   extension. If the starting structure is already
%                   extended with other size, the spins in the original
%                   unit cell are used. Magnetic ordering wavevector k will
%                   be set to zero. To generate structure with non-zero k,
%                   use 'helical' or 'direct' option. (default)
%
%                helical:
%                   generates helical structure, starting structure
%                   is defined by param.S, the normal vector of rotation is
%                   stored in param.n, the ordering wavevector is stored in
%                   param.k. If param.S is complex, it is used as basis
%                   vectors to generate magnetic structure according to the
%                   following formula (param.n value then omitted):
%
%                   M_i(r) = Re(Psi_i)cos(2*pi*km*r)+Im(Psi)sin(2*pi*km*r).
%
%                   param.S has to contain either 1 spin direction or basis
%                   vector, or as many as the number of magnetic atoms in
%                   the crystallographic unit cell. In the first case, the
%                   r position is the atomic position, in the second case r
%                   is the lattice translation vector of the
%                   crystallographic cell where the moment directions are
%                   calculated.
%
%                rotate:
%                   uniform rotation of all magnetic moments with a
%                   param.phi angle around the param.n vector. If
%                   param.phi=0, all moments are rotated so, that the first
%                   moment is parallel to param.n vector in case of
%                   collinear structure or in case of planar structure
%                   param.n defines the normal of the plane of the magnetic
%                   moments.
%
%                func:
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
%
% phi       Angle of rotation of the magnetic moments in rad. Default
%           value 0.
%
% nExt      Number of unit cell to extend the magnetic structure,
%           dimensions are [1 3]. Default value is stored in obj.
%
% k         Magnetic ordering wavevector in r.l.u., dimensions are [1 3].
%           Default value is defined in obj.
%
% n         Normal vector to the spin rotation plane, dimensions are [1 3].
%           Default value [0 0 1].
%
% S         Direct input of the spin values, dimensions are [3 nSpin].
%           Default value is stored in swobj.
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
%
% See also SW, SW.ANNEAL, SW.OPTMAGSTR, GM_SPHERICAL3D, GM_PLANAR.
%
%

inpForm.fname  = {'mode'   'nExt'            'k'           'n'     'S'           'phi' 'epsilon'};
inpForm.defval = {'extend' obj.mag_str.N_ext obj.mag_str.k [0 0 1]  obj.mag_str.S 0     1e-5     };
inpForm.size   = {[1 -1]   [1 3]             [1 3]         [1 3]   [3 -2]        [1 1] [1 1]    };
inpForm.soft   = {false    false             false         false   false         false false    };

inpForm.fname  = [inpForm.fname  {'func'          'x0'  }];
inpForm.defval = [inpForm.defval {@gm_spherical3d []    }];
inpForm.size   = [inpForm.size   {[1 1]           [1 -3]}];
inpForm.soft   = [inpForm.soft   {false           true  }];

param = sw_readparam(inpForm, varargin{:});

nExt     = double(param.nExt);
mAtom    = obj.matom;
nMagAtom = size(mAtom.r,2);
nMagExt  = nMagAtom*prod(nExt);

% Create mAtom.Sext matrix.
mAtom    = sw_extendlattice(nExt, mAtom);

% Magnetic ordering wavevector
k = param.k;

% Axis of rotation, size (1,3)
n = (param.n)/norm(param.n);

% If the magnetic structure is not initialized start with a random one.
if strcmp(param.mode,'extend') && (nMagAtom > size(param.S,2))
    param.mode = 'random';
    warning('sw:genmagstr:WrongInitialStructure','No magnetic structure is defined, random structure is created instead!')
end

switch param.mode
    case 'extend'
        % Extend the unit cell if:
        % -the new number of extended cells does not equal to the number of
        %  cells defined in swobj
        % -the number of spins stored in swobj is not equal to the number
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
        S = bsxfun(@times,S,mAtom.Sext')';
        k = [0 0 0];
    case 'helical'
        S0 = param.S;
        % Magnetic ordering wavevector in the extended unit cell.
        kExt = k.*nExt;
        % Warns about the non sufficient extension of the unit cell.
        if any(abs(kExt-round(kExt))>param.epsilon)
            warning('spinw:sw_magstr:UCExtNonSuff','In the extended unit cell k is still larger than epsilon!');
        end
        % Number of spins in the input.
        nSpin = size(param.S,2);
        
        if (nSpin~= nMagAtom) && (nSpin==1)
            % Single defined spin, use the atomic position.
            r = mAtom.RRext;
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
        if isreal(param.S)
            % Rotate spins for each unit cell.
            for ii = 1:nMagExt
                selS    = S0(:,mod(ii-1,nSpin)+1);
                S(:,ii) = sw_rot(n,phi(ii),selS);
            end
        else
            % Use param.S as complex basis vectors.
            for ii = 1:nMagExt
                selS    = S0(:,mod(ii-1,nSpin)+1);
                S(:,ii) = real(selS)*cos(phi(ii))+imag(selS)*sin(phi(ii));
            end
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
            phi = atan2(norm(cross(S1,nRot)),dot(S1,nRot));
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
        
        [S, k, n] = param.func(S, param.x0);
        
    otherwise
        error('sw:genmagstr:WrongMode','Wrong param.mode value!');
end

mag_str.N_ext = int32(nExt(:))';
mag_str.k     = k(:)';
mag_str.S     = S;
mag_str.n     = n(:)';

obj.mag_str   = mag_str;

end