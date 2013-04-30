function varargout = plot(obj, varargin)
% plots crystal structure, magnetic structure, anisotropy and couplings
%
% [hFigure, handle] = PLOT(obj, 'option1', value1, ...): plot the atoms and
% couplings of obj. Instead of ('option1', value1) pairs, a struct can be
% used.
%
% Options:
%
%   range           Plotting range of the lattice parameters, size: (3,2)
%   plotSpin        Whether to plot magnetic structure, default is true.
%   plotC           Plot couplings. Default is true.
%   plotDM          If non-zero, Dzyaloshinskii-Moriya vectors are plotted
%                   onto the middle of the bond. plotDM value is used to
%                   scale the length of the vector. Standard bond is
%                   plotted if plotDM is zero. Default is zero.
%   plotZeroC       Plot couplings with zero value, default is true.
%   legend          Whether to plot legend, default is true.
%   plotCell        Whether to plot unit cell at the origin, defult is true.
%   plotMultCell    Whether to plot multiple unit cells, default is false.
%   plotAxis        Whether to plot (a,b,c) axis, default is true.
%   colorAxis       Color of the (a,b,c) axis, size is (1,3), default is
%                   black ([0 0 0]).
%   surfRes         Number of points on the surface mesh, default is 30.
%   couplingR       Radius of the cylinder of the couplings, default is
%                   0.05 Angstrom.
%   scaleS          Scaling factor for the lenght of the spins, default is
%                   1.
%   radiusS         Radius of the cylinder of the spins, default is 0.06.
%   headAngleS      Angle of the spin arrow head, default is 15 deg.
%   headLengthS     Length of the spin arrow head, default is 0.5
%                   Angstrom.
%   alphaPlaneS     Transparency of the circle, representing rotation
%                   plane of the moments, default is 0.07.
%   atomLabel       Whether to plot labels for atoms, default is true.
%   colorS          Color of the magnetic moment vectors, default is
%                   green ([0 255 0]).
%   colorCell       Color of the unit cell, default is black ([0 0 0]).
%   lineStyleCell   Line style of the unit cell, default is '--'.
%   dAxis           Distance of coordinate system arrows origin from
%                   crystal origin in Angstrom, default is [0.5;1.5;2.0].
%   dText           Distance between item and its text label, default is
%                   0.2 Angstrom.
%   wSpace          Space between figure window and plot structure,
%                   default is 10.
%   rAtom           The default atom radius, default is 0.3
%                   Angstrom.
%   rAtomData       The source of atomic radius data:
%                    false: use rAtom for all atoms (default)
%                    true:  use radius data from database based on atom
%                           label multiplied by rAtom value.
%   alphaAniso      Transparency of the anisotropy sphere, default is 0.3.
%   coplanar        Decides the limits for plotting the plane of the
%                   magnetic moments.
%                    0    No plane is plotted.
%                    x    Plane of the moment is plotted if the best
%                         fitting plane is better than x.
%                   Default is 0.1.
%   ellAniso        Whether to draw the ellipsoid giving the anisotropy.
%                    0:   No ellipsoid is drawn.
%                    R:   The maximum radius of the ellipsoid.
%                   Default is 0.
%   ellEpsilon      Multiplicator for the small axes of the anisotropy
%                   ellipsoid. Default is 0.1.
%   tooltip         Whether to write tooltip information for graphical
%                   objects. Default is true.
%   fontSize        Font size of all text on the plot in pt, default is 12.
%
% handle            Stores the handles for all graphical object in a
%                   struct.
%
% To access the graphical objects, use the function:
%   handleList = sw_getobject(tagName,fHandle)
% This command returns all objects with tagName string in their 'Tag' field,
% in figure with fHandle.
%
% List of Tag values:
% unitCell         Unit cell.
% arrow            Arrows of the lattice vectors.
% arrowText        Text of the lattice vectors.
% atom_            Spheres of the atoms, for example 'atom_Cr' for atoms
%                  with the label 'Cr'.
% aniso_           Anisotropy ellipsoid with atom names, for example
%                  'aniso_Cr', symbolising the anisotropy matrix.
% atomText_        Label of atoms, for example 'atomText_Cr'.
% coupling_        Cylinder of couplings, for example 'coupling_J1' for the
%                  couplings with label 'J1'.
% spinPlane        Circle showing the plane of the spins.
% legend           Legend.
%

lattice   = obj.lattice;
% The basis vectors in columns.
basisVector = obj.basisvector;

% Handle of figure window to plot.
hFigure = sw_getfighandle('sw_crystal');
if isempty(hFigure)
    hFigure = sw_structfigure;
end

tooltip(hFigure,'<< Click again on any object! >>')

cva = get(gca,'CameraViewAngle');
[az, el] = view();

set(gca,'Position',[0 0 1 1]);
set(gca,'Color','none');
set(gca,'Box','off');
set(gca,'Clipping','Off');
daspect([1 1 1]);
pbaspect([1 1 1]);
axis off
axis vis3d
hold on
material dull;

%% input parameters

inpForm.fname  = {'range'                      'plotSpin' 'colorAxis' 'plotCell' 'plotDM' };
inpForm.defval = {[-0.1 1.1;-0.1 1.1;-0.1 1.1] true       [0 0 0]     true       0        };
inpForm.size   = {[3 2]                        [1 1]      [1 3]       [1 1]      [1 1]    };

inpForm.fname  = [inpForm.fname  {'plotMultCell' 'plotAxis' 'surfRes' 'couplingR' 'scaleS' 'radiusS'}];
inpForm.defval = [inpForm.defval {false          true       30        0.05        1        0.06     }];
inpForm.size   = [inpForm.size   {[1 1]          [1 1]      [1 1]     [1 1]       [1 1]    [1 1]    }];

inpForm.fname  = [inpForm.fname  {'headAngleS' 'headLengthS' 'alphaPlaneS'}];
inpForm.defval = [inpForm.defval {15           0.5           0.07         }];
inpForm.size   = [inpForm.size   {[1 1]        [1 1]         [1 1]        }];

inpForm.fname  = [inpForm.fname  {'atomLabel'  'colorS'  'colorCell' 'plotC' 'legend'}];
inpForm.defval = [inpForm.defval {true         [0 255 0] [0 0 0]    true    true    }];
inpForm.size   = [inpForm.size   {[1 1]        [1 3]     [1 3]       [1 1]   [1 1]   }];

inpForm.fname  = [inpForm.fname  {'lineStyleCell' 'dAxis'       'dText' 'wSpace' 'rAtomData'}];
inpForm.defval = [inpForm.defval {'--'            [0.5;1.5;2.0] 0.2     10       false      }];
inpForm.size   = [inpForm.size   {[1 -1]          [3 1]         [1 1]   [1 1]    [1 1]      }];

inpForm.fname  = [inpForm.fname  {'rAtom' 'alphaAniso' 'coplanar' 'fontSize'}];
inpForm.defval = [inpForm.defval {0.3     0.3          0.1        12        }];
inpForm.size   = [inpForm.size   {[1 1]   [1 1]        [1 1]      [1 1]     }];

inpForm.fname  = [inpForm.fname  {'ellAniso' 'ellEpsilon' 'plotZeroC' 'tooltip'}];
inpForm.defval = [inpForm.defval {0          0.1          true        true     }];
inpForm.size   = [inpForm.size   {[1 1]      [1 1]        [1 1]       [1 1]    }];

param = sw_readparam(inpForm, varargin{:});

%shift of the plot
param.wShift = -sum(basisVector * sum(param.range,2)/2,2);
param.wShift = repmat(param.wShift,[1 2]);
param.wShift = reshape(param.wShift',1,[]);

% Plot spins if spin variables exist and param.plotSpin is true.
param.plotSpin = param.plotSpin && ~isempty(obj.mag_str.S);

%% figure data

uData = get(hFigure,'UserData');

uData.obj   = copy(obj);
uData.param = param;

%% Calculated parameters

% basis vectors for unit cell
v0 = zeros(3,1);
v1 = basisVector(:,1);
v2 = basisVector(:,2);
v3 = basisVector(:,3);

%axis range to show
param.wRange = [param.range(:,1)-param.wSpace param.range(:,2)+param.wSpace];
axis(reshape((basisVector*param.wRange)',1,[])+param.wShift);
zoom(10);

%sphare surface for atoms
[sp.x, sp.y, sp.z] = sphere(param.surfRes);

%% Plot the unit cell.

handle.unitCell = [];

if param.plotCell
    
    path=[v0 v1 v1+v2 v2 v0 v3 v3+v1 v1 v3+v1 v3+v1+v2 v1+v2 v3+v1+v2 v3+v2 v2 v3+v2 v3];
    
    if param.plotMultCell
        for ii = ceil(param.range(1,1)):floor(param.range(1,2)-1)
            for jj = ceil(param.range(2,1)):floor(param.range(2,2)-1)
                for kk = ceil(param.range(3,1)):floor(param.range(3,2)-1)
                    pathd = path+(v1*ii+v2*jj+v3*kk)*ones(1,16);
                    handle.unitCell(end+1) = plot3(pathd(1,:),pathd(2,:),pathd(3,:));
                end
            end
        end
    end
    if isempty(handle.unitCell)
        handle.unitCell = plot3(path(1,:),path(2,:),path(3,:));
    end
    
    set(handle.unitCell,'Color',param.colorCell,'LineStyle',param.lineStyleCell);
    set(handle.unitCell,'Tag','unitCell');
    tooltip(handle.unitCell,['Unit cell \na=' sprintf('%5.3f',lattice.lat_const(1)) ' \nb=' sprintf('%5.3f',lattice.lat_const(2)) ' \nc=' sprintf('%5.3f',lattice.lat_const(3))])
end

%% Plot the coordinate system.
if param.plotAxis
    
    r = -ones(3,1).*param.dAxis + v1*param.range(1,1) + v2*param.range(2,1) + v3*param.range(3,1);
    handle.arrow = zeros(12,1);
    
    [handle.arrow(1:4)]  = sw_arrow(r,(r+v1),param.radiusS,param.headAngleS,param.headLengthS,param.surfRes);
    [handle.arrow(5:8)]  = sw_arrow(r,(r+v2),param.radiusS,param.headAngleS,param.headLengthS,param.surfRes);
    [handle.arrow(9:12)] = sw_arrow(r,(r+v3),param.radiusS,param.headAngleS,param.headLengthS,param.surfRes);
    
    set(handle.arrow,'Facecolor',param.colorAxis);
    set(handle.arrow,'Tag','arrow');
    
    v1n = v1/norm(v1)*param.dText;
    v2n = v2/norm(v2)*param.dText;
    v3n = v3/norm(v3)*param.dText;
    
    handle.arrowText(1) = text(r(1)+v1(1)+v1n(1),r(2)+v1(2)+v1n(2),r(3)+v1(3)+v1n(1),'a','fontSize',param.fontSize);
    handle.arrowText(2) = text(r(1)+v2(1)+v2n(1),r(2)+v2(2)+v2n(2),r(3)+v2(3)+v2n(2),'b','fontSize',param.fontSize);
    handle.arrowText(3) = text(r(1)+v3(1)+v3n(1),r(2)+v3(2)+v3n(2),r(3)+v3(3)+v3n(3),'c','fontSize',param.fontSize);
    
    set(handle.arrowText,'Tag','arrowText');
    tooltip(handle.arrow,'Arrow')
    tooltip(handle.arrowText,'Arrow text')
end

handle.atom        = [];
handle.aniso       = [];
handle.atomText    = [];
handle.spinCircle  = [];
handle.spinCircle2 = [];
handle.spinArrow   = [];
handle.coupling    = [];
handle.DMcoupling  = [];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ATOMS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

atom       = obj.atom;
mAtom      = obj.matom;

nAtom      = size(atom.r,2);
nMagAtom   = size(mAtom.r,2);

atom.color = double(obj.unit_cell.color(:,atom.idx))/255;
atom.name  = obj.unit_cell.label(atom.idx);
atom.rad   = ones(nAtom,1) * param.rAtom; % atomic radius

for ll = 1:nAtom
    % Saves atomic radius data if allowed and exists scaled with param.rAtom.
    if param.rAtomData && (sw_ratom(atom.name{ll}) > 0)
        atom.rad(ll) = sw_ratom(atom.name{ll})*param.rAtom;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% COUPLINGS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

matrix   = obj.matrix;
coupling = obj.coupling;
nExt     = double(obj.mag_str.N_ext);

cIdx = any(coupling.mat_idx,1);
coupling.dl      = coupling.dl(:,cIdx);
coupling.atom1   = coupling.atom1(cIdx);
coupling.atom2   = coupling.atom2(cIdx);
coupling.mat_idx = coupling.mat_idx(:,cIdx);
coupling.idx   = coupling.idx(cIdx);

nCoupling     = size(coupling.idx,2);

% Create the actual indexing and type in coupling
coupling.idxJ = zeros(3,nCoupling);
coupling.type = zeros(3,nCoupling);

for ii = 1: nCoupling
    for jj = 1:3
        idxJ = coupling.mat_idx(jj,ii);
        if any(idxJ)
            % save zero valued couplings if plotZeroC is true
            if param.plotZeroC || any(any(matrix.mat(:,:,idxJ)))
                coupling.idxJ(jj,ii) = idxJ;
                coupling.type(jj,ii) = sw_mattype(matrix.mat(:,:,idxJ));
            end
            
        end
    end
end




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MAGNETIC MOMENTS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if param.coplanar>0
    [nSpin, param.coplanar] = coplanar([[0;0;0] obj.mag_str.S],param.coplanar);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ANISOTROPY ELLIPSOIDS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

single_ion  = obj.single_ion;
aniso.mat   = zeros(3,3,nMagAtom);
aniso.ell   = zeros(3,3,nMagAtom);
aniso.label = cell(1,nMagAtom);

propAniso = true;

if param.ellAniso>0
    % Creates the anisotropy ellipsoids.
    if length(single_ion.aniso) ~= nMagAtom
        error('sw:plot:AnisotropyMissing','sw.single_ion field is not initialized, use: sw.gencoupling!');
    end
    for ll = 1:nMagAtom
        % Selects the anisotropy matrix
        cIdx = single_ion.aniso(ll);
        if any(cIdx)
            aniso.mat(:,:,ll) = matrix.mat(:,:,cIdx);
            
            
            % Calculates the main radiuses of the ellipsoid.
            [V, R] = eig(aniso.mat(:,:,ll));
            if norm(aniso.mat(:,:,ll)-aniso.mat(:,:,ll)') > 1e-8
                propAniso = false;
            end
            % Creates positive definite matrix by adding constant to all
            % eigenvalues.
            R0      = diag(R);
            epsilon = sqrt(sum(R0.^2))*param.ellEpsilon;
            dR      = 1./(R0-min(R0)+epsilon);
            dR0     = max(dR);
            if dR0 ~= 0
                dR = dR/dR0;
            else
                dR = [1 1 1];
            end
            
            R = diag(dR)*param.ellAniso;
            
            aniso.ell(:,:,ll) = V*R;
            aniso.label{ll}   = matrix.label{cIdx};
        end
    end
end

if ~propAniso
    warning('sw:plot:WrongAnisotropy','Anisotropy matrix is asymmetric, cannot plot ellipsoid!');
end


%% Plot all atoms and couplings in selected range.
for ii = floor(param.range(1,1)):floor(param.range(1,2))
    for jj = floor(param.range(2,1)):floor(param.range(2,2))
        for kk = floor(param.range(3,1)):floor(param.range(3,2))
            dCell = [ii;jj;kk];
            % Store the index of magAtom
            llMagAtom = 1;
            for ll = 1:nAtom
                % Index of spin in the extended unit cell
                llSpin = mod(dCell,nExt')'*[1; nExt(1); nExt(1)*nExt(2)]*nMagAtom;
                
                % Atomic position in lattice units.
                rLat = atom.r(:,ll) + dCell;
                
                % Plot if the atom is within the plotting range.
                if all((rLat < param.range(:,2)) & (rLat > param.range(:,1)))
                    rPlot = basisVector*rLat;
                    
                    % Draw anisotropy semi-transparent ellipsoid.
                    if atom.mag(ll) && param.ellAniso>0 && propAniso
                        aIdx = obj.single_ion.aniso(llMagAtom);
                        if any(aIdx) && (norm(obj.matrix.mat(:,:,aIdx))>0)
                            
                            ell.xyz = aniso.ell(:,:,llMagAtom)*[sp.x(:) sp.y(:) sp.z(:)]';
                            
                            ell.x = reshape(ell.xyz(1,:),[1 1]*param.surfRes+1);
                            ell.y = reshape(ell.xyz(2,:),[1 1]*param.surfRes+1);
                            ell.z = reshape(ell.xyz(3,:),[1 1]*param.surfRes+1);
                            
                            sAniso = surf(ell.x+rPlot(1),ell.y+rPlot(2),ell.z+rPlot(3));
                            
                            handle.aniso(atom.idx(ll),end+1) = sAniso;
                            set(sAniso,'LineStyle','none');
                            set(sAniso,'FaceAlpha',param.alphaAniso);
                            set(sAniso,'FaceColor',double(obj.matrix.color(:,aIdx))/255);
                            set(sAniso,'Tag',['aniso_' atom.name{ll}]);
                            
                            tooltip(sAniso,[aniso.label{llMagAtom} ' anisotropy\n' strmat(aniso.mat(:,:,llMagAtom))]);
                        end
                    end
                    
                    % Plot the label of the atom.
                    if (~any(dCell)) && param.atomLabel
                        
                        handle.atomText(atom.idx(ll),end+1) = text(...
                            'Position',             rPlot+param.dText+atom.rad(ll)/sqrt(3),...
                            'String',               sprintf('%s(%d)_%d',atom.name{ll},atom.idx(ll),ll),...
                            'VerticalAlignment',    'bottom',...
                            'Tag',                  ['atomText_' atom.name{ll} ' '],...
                            'fontSize',             param.fontSize);
                        tooltip(handle.atomText(atom.idx(ll),end),[atom.name{ll} ' atom \nUnit cell: \n' sprintf('[%d, %d, %d]',dCell) '\nAtomic position: \n' sprintf('[%6.3f, %6.3f, %6.3f] ',atom.r(:,ll))]);
                    end
                    
                    % Creates the sphere of the atom.
                    aSphere = surf(sp.x*atom.rad(ll)+rPlot(1), sp.y*atom.rad(ll)+rPlot(2), sp.z*atom.rad(ll)+rPlot(3));
                    set(aSphere,'Tag',['atom_' atom.name{ll}]);
                    handle.atom(atom.idx(ll),end+1) = aSphere;
                    set(aSphere,'LineStyle','none');
                    set(aSphere,'FaceColor',atom.color(:,ll));
                    
                    tooltip(aSphere,[atom.name{ll} ' atom \nUnit cell: \n' sprintf('[%d, %d, %d]',dCell) '\nAtomic position: \n' sprintf('[%6.3f, %6.3f, %6.3f] ',atom.r(:,ll))]);
                    
                    % Plots magnetic moments.
                    if param.plotSpin && atom.mag(ll)
                        
                        plotS   = param.scaleS*obj.mag_str.S(:,llMagAtom+llSpin);
                        lengthS = norm(obj.mag_str.S(:,llMagAtom+llSpin));
                        
                        if param.coplanar
                            % TODO
                            % Circle for planar structure and maybe cones :)
                            cPoints = sw_circle(rPlot,nSpin,norm(plotS),param.surfRes);
                            C1 = plot3(cPoints(1,:),cPoints(2,:),cPoints(3,:));
                            handle.spinCircle(atom.idx(ll),end+1) = C1;
                            set(C1,'Color',atom.color(:,ll));
                            set(C1,'Tag','spinPlane');
                            tooltip(C1,['Plane of the ordered moment: \n' sprintf('[%6.3f, %6.3f, %6.3f]',nSpin)]);
                            
                            C2 = sw_circlesurf(rPlot,nSpin,norm(plotS),param.surfRes);
                            handle.spinCircle2(atom.idx(ll),end+1) = C2;
                            set(C2,'LineStyle','none');
                            set(C2,'FaceAlpha',param.alphaPlaneS);
                            set(C2,'FaceColor',atom.color(:,ll));
                            tooltip(C2,['Plane of the ordered moment: \n' sprintf('[%6.3f, %6.3f, %6.3f]',nSpin)]);
                        end
                        
                        halfS = false;
                        intS  = false;
                        if abs(round(lengthS)-lengthS)<0.01
                            lengthS = round(lengthS);
                            halfS   = false;
                            intS    = true;
                        elseif abs(round(lengthS*2)-lengthS*2)<0.02
                            lengthS2 = round(lengthS*2);
                            halfS = true;
                            intS  = true;
                        end
                        
                        
                        hArrow  = sw_arrow(rPlot-plotS,rPlot+plotS,param.radiusS,param.headAngleS,param.headLengthS,param.surfRes);
                        set(hArrow,'Tag','spinArrow');
                        handle.spinArrow(atom.idx(ll),end+(1:4)) = hArrow;
                        set(hArrow,'FaceColor',param.colorS/255);
                        set(hArrow,'LineStyle','none');
                        if intS
                            if halfS
                                tooltip(hArrow,['S=' sprintf('%d',lengthS2) '/2 magnetic moment \nxyz components:  \n' sprintf('[%6.3f, %6.3f, %6.3f]',obj.mag_str.S(:,llMagAtom+llSpin))])
                            else
                                tooltip(hArrow,['S=' sprintf('%d',lengthS) ' magnetic moment \nxyz components:  \n' sprintf('[%6.3f, %6.3f, %6.3f]',obj.mag_str.S(:,llMagAtom+llSpin))])
                            end
                        else
                            tooltip(hArrow,['S=' sprintf('%4.2f',lengthS) ' magnetic moment \nxyz components:  \n' sprintf('[%6.3f, %6.3f, %6.3f]',obj.mag_str.S(:,llMagAtom+llSpin))])
                        end
                    end
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                end
                if atom.mag(ll)
                    llMagAtom = llMagAtom + 1;
                end
            end
        end
    end
end

% Plots cylinder for couplings between amgnetic atoms.

for ii = floor(param.range(1,1)):floor(param.range(1,2))
    for jj = floor(param.range(2,1)):floor(param.range(2,2))
        for kk = floor(param.range(3,1)):floor(param.range(3,2))
            dCell = [ii;jj;kk];
            for ll = 1:nCoupling
                for mm = 1:3
                    if coupling.idxJ(mm,ll)
                        rLat1 = mAtom.r(:,coupling.atom1(ll)) + dCell;
                        rLat2 = mAtom.r(:,coupling.atom2(ll)) + double(coupling.dl(:,ll)) + dCell;
                        
                        if all((rLat1<param.range(:,2))&(rLat1>param.range(:,1))&(rLat2<param.range(:,2))&(rLat2>param.range(:,1)))
                            rPlot1    = basisVector*rLat1;
                            rPlot2    = basisVector*rLat2;
                            cColor    = double(matrix.color(:,coupling.idxJ(mm,ll)))/255;
                            
                            % plot normal cylinder
                            if param.plotC
                                cCylinder = sw_cylinder(rPlot1,rPlot2,param.couplingR,param.surfRes);
                                handle.coupling(coupling.idx(ll),end+1) = cCylinder;
                                set(cCylinder,'FaceColor',cColor);
                                set(cCylinder,'Tag',['coupling_' matrix.label{coupling.idxJ(mm,ll)}]);
                                
                                tooltip(cCylinder,[matrix.label{coupling.idxJ(mm,ll)} ' coupling\nValue:\n' strmat(matrix.mat(:,:,coupling.idxJ(mm,ll)))]);
                            end
                            % plot DM vector
                            if param.plotDM>0 && (coupling.type(mm,ll)==3)
                                rCent = (rPlot1+rPlot2)/2;
                                mDM   = obj.matrix.mat(:,:,coupling.idxJ(mm,ll));
                                vDM   = [mDM(2,3); mDM(3,1); mDM(1,2)]*param.plotDM;
                                hDM   = sw_arrow(rCent,rCent+vDM,param.radiusS,param.headAngleS,param.headLengthS,param.surfRes);
                                handle.DMcoupling(coupling.idx(ll),end+(1:4)) = hDM;
                                set(hDM,'FaceColor',cColor);
                                set(hDM,'Tag',['coupling_' matrix.label{coupling.idxJ(mm,ll)}]);
                                
                                tooltip(hDM,[matrix.label{coupling.idxJ(mm,ll)} ' DM coupling\nValue:\n' strmat(vDM')]);
                                if ~param.plotC
                                    rPlot = [rPlot1 rPlot2];
                                    hLine = line(rPlot(1,:),rPlot(2,:),rPlot(3,:),'LineStyle','--','LineWidth',2,'Color',[0 0 0]);
                                    set(hLine,'Tag',['coupling_' matrix.label{coupling.idxJ(mm,ll)}]);
                                    handle.DMcoupling(coupling.idx(ll),end+1) = hLine;
                                end
                            end
                        end
                    end
                end
            end
        end
    end
end


%% Do the rest.

view([az el]);
hold off

if isfield(uData,'handle')
    oldHandle = uData.handle;
    hName = fieldnames(oldHandle);
    
    for ii = 1:length(hName)
        h0 = reshape(oldHandle.(hName{ii}),1,[]);
        h0(h0 == 0) = [];
        h0(~ishandle(h0)) = [];
        delete(h0);
    end
end

h2 = hgtransform('Parent',uData.h);
hName = fieldnames(handle);

for ii = 1:length(hName)
    h0 = reshape(handle.(hName{ii}),1,[]);
    h0(h0 == 0) = [];
    h0(~ishandle(h0)) = [];
    set(h0,'Parent',h2);
    set(h0,'Clipping','Off');
end

set(gca,'CameraViewAngle',cva);
handle.light = camlight('right');

T = makehgtform('translate',-sum(basisVector * sum(param.range,2)/2,2)');
set(h2,'Matrix',T);

%% Plot the legend.

cAxis = gca;
if (param.legend) && size(matrix.mat,3)>0
    lHeight = length(matrix.label)*20;
    lWidth  = 80;
    sRadius = 8;
    dA      = 2;
    
    lAxes = axes('Units','pixel','Position',[dA dA lWidth+dA lHeight+dA]);
    axis([0 lWidth 0 lHeight]);
    axis off
    handle.legend = rectangle('Position',[1 1 lWidth-2 lHeight-2],'FaceColor','w');
    set(handle.legend,'Tag','legend');
    hold on
    for ii = 1:length(matrix.label)
        handle.lRect(ii) = rectangle('Position',[5 lHeight-20*ii+6 sRadius*2 sRadius]);
        set(handle.lRect(ii),'FaceColor',double(matrix.color(:,ii))/255);
        handle.lText(ii) = text(30,(lHeight-20*ii+10),matrix.label{ii},...
            'fontSize',param.fontSize);
    end
    set(handle.lRect,'Tag','legend');
    set(handle.lText,'Tag','legend');
    hold off
    
    set(handle.legend,'hittest','off');
    set(handle.lRect,'hittest','off');
    set(handle.lText,'hittest','off');
    set(lAxes,'hittest','off')
    
end
% Tooltip text box
if param.tooltip
    axes('Units','normalized','Position',[0.01 0.9 0.1 0.2]);
    axis off
    handle.tooltip = text(0.05,0.45,'',...
        'units','normalized','horizontalalignment','left','tag','tooltip','FontSize',param.fontSize,'VerticalAlignment','top');
    % Avoid object get active for strange zooming effect
    set(handle.tooltip,'hittest','off');
    
end

axes(cAxis); %#ok<MAXES>

%% Stuff...

uData.handle = handle;

set(hFigure,'UserData',uData);

if get(gca,'CameraViewAngle') == 0.6
    camva('auto');
end;


%% Output

switch nargout
    case 1
        varargout{1} = hFigure;
    case 2
        varargout{1} = hFigure;
        varargout{2} = handle;
end

end

function [n, coplanar] = coplanar(points, epsilon)
% [n, coplanar] = COPLANAR(points,{epsilon}) determines whether the points
% are coplanar. They are coplanar if the mean square of the distances of
% the points from the best fitting plane is smaller than epsilon.
% epsilon       Limit for coplanar structures, default is 0.1.
%
% n             Best normal vector, size is (3,1).
% coplanar      Logical variable.
%

if nargin == 1
    epsilon = 0.1;
end

x = points(1,:)';
y = points(2,:)';
z = points(3,:)';

V = [x-mean(x),y-mean(y),z-mean(z)];
A = (1/length(x))*(V')*V;
[U,D] = svd(A);

n = U(:,end);

coplanar = D(end,end) < epsilon;

end

function string = strmat(mat, format)
% strmat(mat, {format}) creates NxM matrix string, plots nice zeros :)
% mat       NxM matrix
% format    Numerical format for sprintf function, 1x2 vector: [p1 p2]. It
%           generates the format %p1.p2f. Default is [6 3], optional.
%

if nargin == 1
    p1  = 6;
    p2  = 3;
else
    p1 = format(1);
    p2 = format(2);
end

string = [];
for ii = 1:size(mat,1)
    string = [string '|']; %#ok<AGROW>
    for jj = 1:size(mat,2)
        string = [string sprintf(['%' num2str(p1) '.' num2str(p2) 'f '],mat(ii,jj))]; %#ok<AGROW>
    end
    string = [string '|\n']; %#ok<AGROW>
end

end

function tooltip(hObject, string)
% Creates tooltip for graphical objects, with string.

set(hObject,'buttondownfcn',['oTemp=findobj(''tag'',''tooltip''); if ~isempty(oTemp) set(oTemp,''string'',sprintf(''' string '''));end']);

end