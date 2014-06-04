function genlattice(obj, varargin)
% generates crystal lattice from given parameters
%
% GENLATTICE(obj, 'option1', value1, 'option2', value2...)
%
% Input:
%
% obj       sw class object.
%
% Options:
%
% angled    Alpha, beta, gamma angles in degree, dimensions are [1 3].
% angle     Alpha, beta, gamma angles in radian, dimensions are [1 3].
% lat_const a, b, c lattice parameters, dimensions are [1 3].
% sym       Space group index, or space group name (string).
%
% Alternatively the lattice parameters can be given directly when the sw
% object is created using: sw(inpStr), where struct contains the fields
% with initial parameters, e.g.:
%   inpStr.lattice.lat_const = [3 3 4];
%
% The sym option points to the appropriate line in the symmetry.dat file,
% where every line defines a space group by its generators. The first 230
% lines contains all crystallographic space groups with standard setting
% as it is in the International Tables of Crystallography. Additional lines
% can be added to the symmetry.dat file using the sw_addsym() function.
% Every line in the symmetry.dat file can be referenced by either its line
% index or by its label (string).
%
% If the sym option is 0, no symmetry will be used. The sw.gencoupling()
% function will determine the equivalent bonds based on bond length.
%
% Output:
%
% The obj.lattice field will be changed based on the input, the lattice
% constants stored directly and the optional space group string is
% converted to the integer type index.
%
% Example:
%
% ...
% crystal.genlattice('lat_const',[3 3 4],'angled',[90 90 120],'sym','P 6')
% crystal.genlattice('lat_const',[3 3 4],'angled',[90 90 120],'sym',168)
%
% The two lines are equivalent, both will create hexagonal lattice, with 
% 'P 6' space group.
%
% See also SW, SW_ADDSYM, SW_GENSYM, SW.GENCOUPLING.
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