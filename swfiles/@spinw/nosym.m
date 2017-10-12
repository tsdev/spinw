function nosym(obj, varargin)
% reduces symmetry to P0
% 
% ### Syntax
% 
% `nosym(obj)`
% 
% ### Description
% 
% `nosym(obj)` reduces the crystal symmetry to $P0$ but keeps all symmetry
% generated atoms, that become all symmetry inequivalent. The function can
% be used to test different types of symmetry breaking terms in the spin
% Hamiltonian.
% 
% ### Examples
% 
% The example generates an FCC cell using explicit translations. After
% applying the `spinw.nosym` function, the `cryst.unit_cell.r` contains the
% four generated atomic positions, that are not symmetry equivalent any
% more.
%
% ```
% >>symOp = 'x+1/2,y+1/2,z;x+1/2,y,z+1/2;x,y+1/2,z+1/2'
% >>cryst = spinw
% >>cryst.genlattice('lat_const',[8 8 8],'sym',symOp,'label','FCC')
% >>cryst.addatom('r',[0 0 0],'label','Atom1')
% >>cryst.unit_cell.r>>
% >>cryst.nosym
% >>cryst.unit_cell.r>>
% ```
% 
% ### Input Arguments
% 
% `obj`
% : [spinw] object.
% 
% ### Output Arguments
% 
% The `obj` input will have `obj.lattice.sym` field equal to zero and the
% [obj.unit_cell] field will contain all the generated atomic positions.
% 
% ### See Also
% 
% [spinw] \| [spinw.newcell]
%

% TODO keep bonds when removing symmetry

obj.newcell('bvect',{[1 0 0] [0 1 0] [0 0 1]});
obj.lattice.sym = zeros(3,4,0);
obj.sym = false;

end