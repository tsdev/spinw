function genlattice(obj, varargin)
% generates lattice from given parameters
%
% GENLATTICE(obj, 'option1', value1, 'option2', value2...)
%
% Options:
%
%   angled      Alpha, beta, gamma angles in degree, dimensions are [1 3].
%   angle       Alpha, beta, gamma angles in radian.
%   lat_const   Lattice parameters, dimensions are [1 3].
%   sym         Crystal symmetry group index, see symmetry.dat file in the
%               sw folder. To include new symmetry operators, append a new
%               line to symmetry.dat using sw_addsym.
%
% The lattice parameters can be given directly when the sw object is
% created using: sw(param), where struct contains the fields with initial
% parameters, e.g.:
%   param.lattice.lat_const = [3 3 4];
%
% Example:
%
% genlattice(crystal,'lat_const',[3 3 4],'angled',[90 90 120],'sym','P 6');
% This line will create hexagonal lattice, with 'P 6' space group.
%
% See also SW, SW_ADDSYM, SW_GENSYM.
%

inpForm.fname  = {'angle'           'lat_const'           'sym'           'angled'};
inpForm.defval = {obj.lattice.angle obj.lattice.lat_const obj.lattice.sym [0 0 0] };
inpForm.size   = {[1 3]             [1 3]                 [1 -1]          [1 3]   };

param = sw_readparam(inpForm, varargin{:});

if all(param.angled)
    param.angle = param.angled*pi/180;
end

obj.lattice.angle     = param.angle;
obj.lattice.lat_const = param.lat_const;

if numel(param.sym) > 1
    % get the number of the symmetry
    [~,~,~,~,param.sym] = sw_gensym(param.sym);
end

obj.lattice.sym       = int32(param.sym);

end