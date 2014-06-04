function nosym(obj, varargin)
% removes the space group symmetry
%
% NOSYM(obj)
%
% The function reduces the crystal symmetry to  0 and but keeps all atomic
% positions, that are now all symmetry inequivalent.
%
% Input:
%
% obj       sw class object.
%
% Output:
%
% The obj input will have obj.lattice.sym field equal to zero and the
% obj.unit_cell field will contain all the generated atomic positions.
%
% Example:
%
% sw_addsym('x+1/2,y+1/2,z;x+1/2,y,z+1/2;x,y+1/2,z+1/2','FCC');
% cryst = sw;
% cryst.genlattice('lat_const',[8 8 8],'sym','FCC')
% cryst.addatom('r',[0 0 0],'label','Atom1')
% cryst.nosym
%
% The example creates four equivalent atomic positions, after the nosym()
% command, the cryst.unit_cell.r contains the four generated positions,
% that are not symmetry equivalent any more.
%
% See also SW, SW.NEWCELL.
%

obj.newcell({[1 0 0] [0 1 0] [0 0 1]});
obj.lattice.sym = int32(0);
obj.sym = false;

end