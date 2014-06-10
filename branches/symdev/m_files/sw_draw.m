function polyDat = sw_draw(varargin)
% polyDat = SW_DRAW(Option1, Value1, ...) extra plots onto the crystal
% structure plot
%
% SW_DRAW plots extra things (bonds, polyhedra) onto the crystal structure
% plot produced previously by sw.plot() function.
%
% Options:
%
% mode      Selects what to plot:
%           'poly'  draws polyhedra around the center atoms.
%           'bond'  draws bonds around the center atoms.
%           Default is 'poly'.
% cAtom     Indices of atom types for the center atom
%           Default is 1. Can be also a string or strings in a cell that
%           identify atoms by their labels.
% pAtom     Indices of atom types of the surrounding atoms
%           Default is 2. Can be also string or cell of strings.
% range     Plot range in reciprocal lattice units, dimensions are [3 2].
%           Default is the plotting range of the figure.
% limit     Can be a single number: gives the number of neighbours or
%           vector: [min max] distance range from neighbours. Default is 6
%           to plot octahedra.
% edge      Whether to paint the edge of the surfaces to the color of
%           the central atom (true) or keep it black (false). Default is
%           true.
% alpha     Transparency of the plotted surfaces. Default is 1 for
%           non-transparecy.
% cBond     Color of different bonds. Default is 'auto', when they are set
%           to the color of the center atom. [R G B] will fix the color of
%           all bonds to a uniform one.
% rBond     Radius of the cylinder of the bonds, default is 0.15 Angstrom.
% surfRes   Number of points on the surface mesh, default is 30.
%           
%
%
% Output:
%
% polyDat is structure type, with the following fields:
%
% surf          Vector, contains the handle of all plotted surfaces.
% index         Index of the center atom.
% pos           Positions of the surrounding atoms relative to the center
%               atom in Angstrom units.
% center        Positions of the center atoms in the crystal, in Angstrom
%               units.
%
% For example:
%
% plot(cryst);
% sw_draw('mode','bond','cAtom','Cr','pAtom','O','range',6)
%
% It will plot 6 bonds between every chromium atom and the 6 closes oxygen
% atoms.
%
% See also SW, SW_T2G, SW_ORBITAL, SW_ADDOBJECT.
%

if nargin == 0
    help sw_draw;
    return
end

hFigure   = sw_getfighandle('sw_crystal');
if isempty(hFigure)
    error('sw:sw_draw:NoFigure','No active crystal structure figure!');
end

obj         = getappdata(hFigure,'obj');
unit_cell   = obj.unit_cell;
basisVector = obj.basisvector;
param = getappdata(hFigure,'param');

inpForm.fname  = {'mode' 'cAtom' 'pAtom' 'range'     'limit' 'edge' 'alpha' 'rBond' 'surfRes'};
inpForm.defval = {'poly' 1       2       param.range 6       true   1       0.15    30       };
inpForm.size   = {[1 -1] [1 -2]  [1 -3]  [3 2]       [1 -4]  [1 1]  [1 1]   [1 1]   [1 1]    };

param  = sw_readparam(inpForm, varargin{:});

       
% identify atoms by their labels
if ischar(param.cAtom)
    param.cAtom = {param.cAtom};
end
if ischar(param.pAtom)
    param.pAtom = {param.pAtom};
end

if iscell(param.cAtom)
    cAtom = [];
    for ii = 1:numel(param.cAtom)
        cAtom = [cAtom find(cellfun(@isempty,strfind(obj.unit_cell.label,param.cAtom{ii}))==false)]; %#ok<AGROW>
    end
    param.cAtom = cAtom;
end

if iscell(param.pAtom)
    pAtom = [];
    for ii = 1:numel(param.pAtom)
        pAtom = [pAtom find(cellfun(@isempty,strfind(obj.unit_cell.label,param.pAtom{ii}))==false)]; %#ok<AGROW>
    end
    param.pAtom = pAtom;
end

limit = param.limit;

atom = obj.atom;

atom1.r0     = atom.r(:,ismember(atom.idx,param.cAtom));
atom2.r0     = atom.r(:,ismember(atom.idx,param.pAtom));
atom1.index0 = atom.idx(ismember(atom.idx,param.cAtom));
atom2.index0 = atom.idx(ismember(atom.idx,param.pAtom));
atom1.color0 = double(unit_cell.color(:,atom1.index0))/255;
atom2.color0 = double(unit_cell.color(:,atom2.index0))/255;

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

if isempty(atom1.r)
    fprintf0(obj.fileid,'Nothing to plot in the plotting range!\n');
    return
end

nPol  = 9*size(atom2.r,2);
dist  = zeros(1,nPol);
pos   = zeros(3,nPol);

switch param.mode
    case 'poly'
        polyDat.surf   = zeros(1,size(atom1.r,2));
    case 'bond'
        polyDat.surf   = [];
    otherwise
        error('sw_draw:WrongInput','Wrong ''mode'' option, check ''help sw_draw''!');
end

polyDat.index  = zeros(1,size(atom1.r,2));
polyDat.pos    = cell(1,size(atom1.r,2));
polyDat.center = basisVector*atom1.r;

hold on

for ii=1:size(atom1.r,2)
    cpos = basisVector*atom1.r(:,ii);
    index=1;
    for jj = (floor(atom1.r(1,ii))+(-1:1))
        for kk = (floor(atom1.r(2,ii))+(-1:1))
            for ll = (floor(atom1.r(3,ii))+(-1:1))
                
                tr = [jj;kk;ll];
                for mm=1:size(atom2.r0,2)
                    %dist(index)  = norm(basisVector*(atom2.r0(:,mm)-atom1.r(:,ii)+tr));
                    dist(index)  = norm(basisVector*(atom2.r0(:,mm)+tr)-cpos);
                    pos(:,index) = basisVector*(atom2.r0(:,mm)+tr);
                    color2(:,index) = atom2.color0(:,mm);
                    index        = index+1;
                end
            end
        end
    end
    [~, newindex] = sortrows(dist',1);
    dist = dist(newindex);
    pos  = pos(:,newindex);
    color2  = color2(:,newindex);
    if length(limit) == 1
        posP = pos(:,1:limit);
        colorP = color2(:,1:limit);
    else
        posP = pos(:,(dist>=limit(1))&(dist<=limit(2)));
        colorP = color2(:,(dist>=limit(1))&(dist<=limit(2)));
    end
    
    switch param.mode
        case 'poly'
            h1 = polyhedron2(posP);
            set(h1,'Tag',sprintf('poly_%d',ii));
            polyDat.pos{ii} = bsxfun(@minus,posP,polyDat.center(:,ii));
            
            set(h1,'FaceColor',atom1.color(:,ii));
            set(h1,'EdgeColor',atom1.color(:,ii));
            set(h1,'FaceAlpha',param.alpha);
            polyDat.surf(ii)  = h1;
            polyDat.index(ii) = atom1.index(ii);
        case 'bond'
            for jj = 1:size(posP,2)
                cCylinder1 = sw_cylinder(cpos,(posP(:,jj)+cpos)/2,param.rBond,param.surfRes,0);
                cCylinder2 = sw_cylinder((posP(:,jj)+cpos)/2,posP(:,jj),param.rBond,param.surfRes,0);
                if ~isempty(cCylinder1)
                    polyDat.surf(end+1:end+(numel(cCylinder1))) = cCylinder1;
                end
                if ~isempty(cCylinder2)
                    polyDat.surf(end+1:end+(numel(cCylinder2))) = cCylinder2;
                end
                set(cCylinder1,'FaceColor',atom1.color(:,ii));
                set(cCylinder1,'EdgeColor',atom1.color(:,ii));
                set(cCylinder2,'FaceColor',colorP(:,jj));
                set(cCylinder2,'EdgeColor',colorP(:,jj));
                polyDat.index(jj,ii) = atom1.index(ii);
                
            end
    end
end

handle.surf = polyDat.surf;
sw_addobject(hFigure,handle);

if ~param.edge
    set(polyDat.surf,'EdgeColor',[0 0 0]);
end

hold off

end

function [handle] = polyhedron2(vertex)
% vertex: 3x3 array, column vectors contain coordinates
% polyhedron from triangular faces
% Euler: face=2*vertex-4
%

try
    % 3D polyhedra
    trids  = convhulln(vertex');
    handle = patch('faces',trids,'vertices',vertex');
catch %#ok<CTCH>
    % flat polygons
    n1 = -cross(vertex(:,2)-vertex(:,1),vertex(:,3)-vertex(:,1));
    n1 = n1/norm(n1);
    n2  = [0; n1(3); -n1(2)];
    n2(3,~any(abs(n2)>1e-10)) = 1;
    n2  = n2/norm(n2);
    % n3 = n1 x n2
    n3  = cross(n1,n2);
    
    T  = [n1 n2 n3];
    vT = T'*vertex;
    vT = vT(2:3,:)';
    % sort the points around the polygon
    pIndex = convexHull(DelaunayTri(vT(:,1),vT(:,2)));
    
    % reorder the edge points
    vertex = vertex(:,pIndex);
    handle = fill3(vertex(1,:),vertex(2,:),vertex(3,:),[0 1 0]);
    
end

end