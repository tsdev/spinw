function genlattice(obj, varargin)
% generates lattice from given parameters
%
% GENLATTICE(obj, 'option1', value1, 'option2', value2...)
%
% Options:
%
%   angle       Alpha, beta, gamma angles, dimensions are [1 3].
%   lat_const   Lattice parameters, dimensions are [1 3].
%   sym         Crystal symmetry group index, see symmetry.dat file in the
%               sw folder. To include new symmetry operators, append a new
%               line to symmetry.dat.
%
% The lattice parameters can be given directly when the sw object is
% created using: sw(param), where struct contains the fields with initial
% parameters, e.g.:
%   param.lattice.lat_const = [3 3 4];
%
% See also SW.
%

inpForm.fname  = {'angle'      'lat_const' 'sym'};
inpForm.defval = {[1 1 1]*pi/2 [3 3 3]     1    };
inpForm.size   = {[1 3]        [1 3]       [1 1]};

param = sw_readparam(inpForm, varargin{:});

obj.lattice.angle     = param.angle;
obj.lattice.lat_const = param.lat_const;
obj.lattice.sym       = int32(param.sym);

end