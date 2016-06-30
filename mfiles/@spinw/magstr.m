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
%nCell = prod(nExt0);
% number of magnetic atoms in the magnetic cell
%nMagExt = size(obj.mag_str.F,2);
% number of magnetic atom in the unit cell
nAtom   = size(obj.matom.r,2);
% number of k-vectors in the magnetic unit cell
nK      = size(obj.mag_str.k,2);
% if the new supercell not equal to the old, tile up the magnetic moments
if any(nExtNew~=nExt0)
    M0 = repmat(reshape(obj.mag_str.F,[3 nAtom nExt0 nK]),[1 1 ceil(nExtNew./nExt0) 1]);
    
    if any(mod(nExtNew./nExt0,1))
        % remove additional moments if nExtNew is not integer multiples of nExt0
        M0 = M0(:,:,1:nExtNew(1),1:nExtNew(2),1:nExtNew(3),:);
    end
else
    M0 = reshape(obj.mag_str.F,[3 nAtom nExt0 nK]);
end

% create the cell indices for all magnetic atoms in the new supercell
nExt1 = nExtNew-1;
[cIdx{1:3}] = ndgrid(0:nExt1(1),0:nExt1(2),0:nExt1(3));
% dimensions: nExtNew(1) x nExtNew(2) x nExtNew(3) x 3
cIdx        = cat(4,cIdx{:});

% calculate the translation vectors that generate the rotations in the new supercell
tIdx = floor(bsxfun(@rdivide,cIdx,permute(nExt0,[1 3 4 2])));
% propagation vector in the original supercell
kExt0 = bsxfun(@times,obj.mag_str.k,nExt0');
% calculate the phases that generate the rotations for the new supercell
phi  = sum(bsxfun(@times,tIdx,permute(kExt0,[3 4 5 1 2])),4);
% complex phase factors
M = real(bsxfun(@times,M0,exp(2*pi*1i*permute(phi,[4 6 1:3 5]))));
% sum up the wave vectors and reshape to standard dimensions
magOut.S = reshape(sum(M,6),3,[]);
% keep only the first non-zero wave vector
kInc = find(sum(mod(kExt0,1) == 0,1)<3);
if ~isempty(kInc)
    magOut.k = obj.mag_str.k(:,kInc(1))';
    n = cross(real(obj.mag_str.F(:,1,kInc(1))),imag(obj.mag_str.F(:,1,kInc(1))));
    % normalize n-vector
    if norm(n) == 0
        magOut.n = [0 0 1];
    else
        magOut.n = n/norm(n);
    end
else
    magOut.k = [0 0 0];
    magOut.n = [0 0 1];
end


% check whether the above calculation gives an exact magnetic structure
nUn = sw_uniquetol(reshape(cross(real(obj.mag_str.F(:,:,kInc)),imag(obj.mag_str.F(:,:,kInc))),3,[]));

% TODO check whether the original structure is not sinusoidal
%

if numel(kInc)>1 || size(nUn,2)>1
    magOut.exact = false;
else
    magOut.exact = true;
end

magOut.N_ext = nExtNew;

end