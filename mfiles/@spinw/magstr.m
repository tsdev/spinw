function magOut = magstr(obj, varargin)
% generates magnetic structure for the rotating frame
%
% magOut = GENMAGSTR(obj, 'Option1', 'Value1', ...)
%
%

nExt0 = double(obj.mag_str.nExt);

inpForm.fname  = {'nExt'};
inpForm.defval = {nExt0 };
inpForm.size   = {[1 3] };
inpForm.soft   = {false };

param = sw_readparam(inpForm, varargin{:});

% size of the supercell where the magnetic structure is generated
nExtNew = param.nExt;

% number of cells in the magnetic supercell
nCell = prod(nExt0);
% number of magnetic atoms in the magnetic cell
nMagExt = size(obj.mag_str.S,2);
% number of magnetic atom in the unit cell
nAtom   = size(obj.matom.r,2);
% number of k-vectors in the magnetic unit cell
nK      = size(obj.mag_str.k,2);
% if the new supercell not equal to the old, tile up the magnetic moments
if any(nExtNew~=nExt0)
    M0 = repmat(reshape(obj.mag_str.S,[3 nAtom*nExt0 nK]),[1 ceil(nExtNew./nExt0) 1]);
    % remove additional moments if nExtNew is not integer multiples of nExt0
    if 
end


% create the cell indices for all magnetic atoms in the original supercell
nExt1 = nExt0-1;
[cIdx{1:3}] = ndgrid(0:nExt1(1),0:nExt1(2),0:nExt1(3));
% dimensions: nExt(1) x nExt(2) x nExt(3) x 3
cIdx        = cat(4,cIdx{:});

% additional phase for each unit cell within the magnetic supercell
phi = zeros(nK,nCell);
for ii = 1:nK
    phi(ii,:) = reshape(sum(bsxfun(@times,2*pi*obj.mag_str.k(:,ii),permute(cIdx,[4 1:3])),1),[1 nCell]);
end


% Warns about the non sufficient extension of the unit cell.
% we substitute random values for symbolic km
skExt = sw_sub1(kExt,'rand');
if any(abs(skExt-round(skExt))>param.epsilon) && prod(nExt) > 1
    warning('sw:genmagstr:UCExtNonSuff','In the extended unit cell k is still larger than epsilon!');
end


% First crystallographic unit cell defined, use only unit cell
% position.
r = bsxfun(@rdivide,floor(bsxfun(@times,mAtom.RRext,nExt')),nExt');


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


end