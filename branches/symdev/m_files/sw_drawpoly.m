function polyDat = sw_drawpoly(varargin)
% polyDat = SW_DRAWPOLY(Option1, Value1, ...) draws polyhedra around atoms.
%
% Options:
%
% cAtom     Indices of atom types for the center atom
%           Default is 1.
% pAtom     Indices of atom types of the surrounding atoms
%           Default is 2.
% range     Plot range in reciprocal lattice units, dimensions are [3 2].
%           Default is the plotting range of the figure.
% limits    Can be a single number: gives the number of neighbours or
%           vector: [min max] distance range from neighbours.
%           Default is 6 to plot octahedra.
% edge      Whether to paint the edge of the polyhedra to the color of
%           the central atom (true) or keep it black (false).
%           Default is true.
% alpha     Transparency of the polyhedra surface.
%           Default is 0.5.
%
% Output:
% polyDat is structure type, with the following fields:
% polyhedron    Vector, contains the handle of all polyhedron
%               surface.
% index         Index of the center atom.
% pos           Positions of the surrounding atoms relatice to the
%               center atom in Angstrom units.
% center        Positions of the center atoms in the crystal, in
%               Angstrom units.
%
% See also SW, SW_T2G, SW_ORBITAL, SW_ADDOBJECT.
%

if nargin == 0
    help sw_drawpoly;
    return
end

hFigure   = sw_getfighandle('sw_crystal');
if isempty(hFigure)
    error('sw:sw_drawpoly:NoFigure','No active crystal structure figure!');
end

obj         = getappdata(hFigure,'obj');
unit_cell   = obj.unit_cell;
basisVector = obj.basisvector;
param = getappdata(hFigure,'param');

inpForm.fname  = {'cAtom' 'pAtom' 'range'              'limits' 'edge' 'alpha'};
inpForm.defval = {1       2       param.range 6        true   0.5    };
inpForm.size   = {[1 -1]  [1 -2]  [3 2]                [1 -3]   [1 1]  [1 1]  };

param  = sw_readparam(inpForm, varargin{:});
limits = param.limits;

atom = obj.atom;

atom1.r0     = atom.r(:,ismember(atom.idx,param.cAtom));
atom2.r0     = atom.r(:,ismember(atom.idx,param.pAtom));
atom1.index0 = atom.idx(ismember(atom.idx,param.cAtom));
atom1.color0 = double(unit_cell.color(:,atom1.index0))/255;

atom1.r     = [];
atom1.color = [];
atom1.index = [];
atom2.r     = [];

fRange = floor(param.range);

for ii = 1:size(atom1.r0,2)
    for jj = fRange(1,1):fRange(1,2)
        for kk = fRange(2,1):fRange(2,2)
            for ll = fRange(3,1):fRange(3,2)
                
                tr   = [jj;kk;ll];
                rrlu = atom1.r0(:,ii)+tr;
                if ~any((rrlu>=param.range(:,2))|(rrlu<=param.range(:,1)))
                    atom1.r(:,end+1)     = rrlu;
                    atom1.color(:,end+1) = atom1.color0(:,ii);
                    atom1.index(end+1)   = atom1.index0(:,ii);
                end
                
            end
        end
    end
end

nPol  = 9*size(atom2.r,2);
dist  = zeros(1,nPol);
pos   = zeros(3,nPol);

polyDat.polyhedron = zeros(1,size(atom1.r,2));
polyDat.index      = zeros(1,size(atom1.r,2));
polyDat.pos        = cell(1,size(atom1.r,2));
polyDat.center     = basisVector*atom1.r;

hold on

for ii=1:size(atom1.r,2)
    index=1;
    for jj = (floor(atom1.r(1,ii))+(-1:1))
        for kk = (floor(atom1.r(2,ii))+(-1:1))
            for ll = (floor(atom1.r(3,ii))+(-1:1))
                
                tr = [jj;kk;ll];
                for mm=1:size(atom2.r0,2)
                    dist(index)  = norm(basisVector*(atom2.r0(:,mm)-atom1.r(:,ii)+tr));
                    pos(:,index) = basisVector*(atom2.r0(:,mm)+tr);
                    index        = index+1;
                end
            end
        end
    end
    [~, newindex] = sortrows(dist',1);
    dist = dist(newindex);
    pos  = pos(:,newindex);
    if length(limits) == 1
        posP = pos(:,1:limits);
    else
        posP = pos(:,(dist>=limits(1))&(dist<=limits(2)));
    end
    h1 = polyhedron2(posP);
    set(h1,'Tag',sprintf('poly_%d',ii));
    polyDat.pos{ii} = bsxfun(@minus,posP,polyDat.center(:,ii));
    
    set(h1,'FaceColor',atom1.color(:,ii));
    set(h1,'EdgeColor',atom1.color(:,ii));
    set(h1,'FaceAlpha',param.alpha);
    polyDat.polyhedron(ii) = h1;
    polyDat.index(ii)      = atom1.index(ii);
end

handle.polyhedron = polyDat.polyhedron;
sw_addobject(hFigure,handle);

if ~param.edge
    set(polyDat.polyhedron,'EdgeColor',[0 0 0]);
end

hold off

end

function [handle] = polyhedron2(vertex)
%vertex: 3x3 array, column vectors contain coordinates
%polyhedron from triangular faces
%Euler: face=2*vertex-4
%

trids  = convhulln(vertex');
handle = patch('faces',trids,'vertices',vertex');

end