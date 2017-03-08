function sym = symmetry(obj)
% true if space group is used to generate couplings
%
% sym = SYMMETRY(obj)
%
% If true, equivalent couplings are generated based on the
% crystal space group and all matrices (interaction, anisotropy
% and g-tensor) are transformed according to the symmetry
% operators. If false, equivalent couplings are generated based
% on bond length, equivalent matrices won't be transformed
% (all identical).
%
% To change it use spinw.gencoupling with the forceNoSym option.
% To remove all symmetry operators use spinw.nosym.
%
% See also SPINW, SPINW.NOSYM, SPINW.GENCOUPLING.
%

sym = obj.sym;

end