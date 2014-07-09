function varargout = plot(obj, varargin)
% plots crystal structure, magnetic structure, anisotropy and couplings
%
% [hFigure, handle] = PLOT(obj, 'option1', value1, ...)
%
% The function plots the atoms and couplings stored in obj onto a Matlab
% figure. Optionally a script output can be produced that creates a jMol
% plot of the crystal structure.
%
% Input:
%
% obj               sw class object.
%
% Options:
%
%   range           Plotting range of the lattice parameters in lattice
%                   units, dimensions are [3 2]. For example to plot the
%                   first unit cell, use: [0 1;0 1;0 1]. Also the number
%                   unit cells can be given along the a, b and c
%                   directions: [2 1 2], that is equivalent to
%                   [0 2;0 1;0 2]. Default is the single unit cell.
%   format          Output format string with the following values:
%                       'plot'  Plot onto Matlab figure. (Default)
%                       'jmol'  Return Jmol script file as a script,
%                               (.spt file).
%
%   Axis & labels =========================================================
%
%   legend          Whether to plot legend, default is true.
%   pCell           Whether to plot unit cell at the origin, defult is
%                   true.
%   pMultCell       Whether to plot multiple unit cells, default is false.
%   pAxis           Whether to plot (a,b,c) axis, default is true.
%   labelAtom       Whether to plot labels for atoms, default is true.
%   tooltip         Whether to write tooltip information for graphical
%                   objects. Default is true.
%   cCell           Color of the unit cell, default is black ([0 0 0]).
%   cAxis           Color of the (a,b,c) axis, [R G B], default is
%                   black ([0 0 0]).
%   lineStyleCell   Line style of the unit cell, default is '--'.
%   lineWidthCell   Line width of the unit cell, default is 1.
%   dAxis           Distance of coordinate system arrows origin from
%                   crystal origin in Angstrom, default is [0.5;1.5;2.0].
%   dText           Distance between item and its text label, default is
%                   0.2 Angstrom.
%   fontSize        Font size of all text on the plot in pt, default is 12.
%   rAxis           Radius of the cylinder of the axis vectors, default is
%                   0.06.
%   angHeadAxis     Angle of the axis arrow head, default is 15 deg.
%   lHeadAxis       Length of the axis arrow head, default is 0.5
%                   Angstrom.
%
%   Atom ==================================================================
%
%   rAtom           The default atom radius, default is 0.3
%                   Angstrom. If zero, no atoms are drawn.
%   pNonMagAtom     Whether to plot non-magnetic atoms, default is true.
%   rAtomData       The source of atomic radius data:
%                    false: use rAtom for all atoms (default)
%                    true:  use radius data from database based on atom
%                           label multiplied by rAtom value.
%   cAtom           Color of the atoms. Default is 'auto', when all atom
%                   gets the color stored in obj.unit_cell. To fix the
%                   color of all atoms to a single color, give [R G B]
%                   vector as option, containing the RGB code of the color.
%
%   Magnetic structure ====================================================
%
%   pSpin           Whether to plot magnetic structure, default is true.
%   cSpin           Color of the magnetic moment vectors, default is 'auto'
%                   when each spin vector has the color of the magnetic
%                   atom it belongs to.
%   sSpin           Scaling factor for the lenght of the magnetic moments,
%                   default is 1.
%   rSpin           Radius of the cylinder of the spins, default is 0.06.
%   angHeadSpin     Angle of the spin arrow head, default is 15 deg.
%   lHeadSpin       Length of the spin arrow head, default is 0.5
%                   Angstrom.
%   aPlaneSpin      Transparency (alpha value) of the circle, representing
%                   rotation plane of the moments, default is 0.07.
%   coplanar        Decides the limits for plotting the plane of the
%                   magnetic moments.
%                    0    No plane is plotted.
%                    x    Plane of the moment is plotted if the best
%                         fitting plane is better than x.
%                   Default value is 0.
%   centerS         If true, the spin vector is centered on the atom, if
%                   false the beggining of the spin vector is on the atom.
%                   Default is true.
%
%   Anisotropy and g-tensor ===============================================
%
%   sEll            Select, whether to plot ellipsoid for anisotropy or
%                   g-tensor:
%                    0:     no ellipsoid is drawn
%                    1:     anisotropy (default)
%                    2:     g-tensor.
%   rEll            Whether to draw the ellipsoid giving the anisotropy or
%                   g-tensor:
%                    0:   No ellipsoid is drawn.
%                    R:   The maximum radius of the ellipsoid.
%                   Default value is 1.
%   aEll            Transparency (alpha value) of the ellipsoid, default is
%                   0.3.
%   eEll            Minimum radius of the principal axes of the ellipse.
%   lwEll           Line width of the main circles surrounding the
%                   ellipsoid, if zero no circles are drawn. Default is 1.
%
%
%   Couplings =============================================================
%
%   pCoupling       Plot couplings. Default is true.
%   pZeroCoupling   Plot couplings with zero value, default is true.
%   pDM             If non-zero, Dzyaloshinskii-Moriya vectors are plotted
%                   onto the middle of the bond. pDM value is used to
%                   scale the length of the vector. Standard bond is
%                   plotted if pDM is zero. Default is 1.
%   cCoupling       Color of different couplings. Default is 'auto', when
%                   they are set to the color stored in the corresponding
%                   obj.matrix. [R G B] will fix the color of all couplings
%                   to a uniform one.
%   cDM             Color of different DM vectors. Default is 'auto', when
%                   they are set to the color stored in the corresponding
%                   obj.matrix. [R G B] will fix the color of all DM
%                   vectors to a uniform one.
%   rCoupling       Radius of the cylinder of the couplings, default is
%                   0.05 Angstrom.
%
%   Other options =========================================================
%
%   surfRes         Number of points on the surface mesh, default is 30.
%   wSpace          Space between figure window and plot structure,
%                   default is 10.
%   hg              Whether to use hgtransform (nice rotation with the
%                   mouse) or default Matlab rotation of 3D objects.
%                   Default is true.
%   hFigure         Handle of the figure to plot. In not given, the figure
%                   window handle is determined automatically.
%   zoom            Number of steps to zoom in/out on the figure. Positive
%                   number gives zoom in, negative numbers give zoom out.
%                   One step equal to the (+/-) button press on the figure
%                   toolbar. Default is zero.
%
% Output:
%
% For plotting:
% hFigure           Handle of the plot figure.
% handle            Stores the handles for all graphical object in a
%                   struct.
%
% For 'jmol' script file output:
% strOut            A single string is returned, that contains a Jmol
%                   script that reproduce the figure as close as possible.
%
%
% The labels on the atom has the following meaning:
%
% AtomLabel(idx1)_idx2:
%   idx1    The index of the atom in obj.unit_cell
%   idx2    The index of the atom in obj.atom
%
% To access the graphical objects, use the function:
%   handleList = sw_getobject(tagName,fHandle);
% This command returns all objects with tagName string in their 'Tag'
% field, in figure with fHandle.
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
% Example:
%
% sq = sw_model('squareAF',1);
% plot(sq);
%
% The example plots the structure, couplings and magnetic structure of the
% square lattice antiferromagnet.
%
% See also SW_DRAW, SW_GETOBJECT, SW_ADDOBJECT.
%

%% input parameters

range0 = [0 1;0 1;0 1];
inpForm.fname  = {'range' 'pSpin'    'cAxis'    'pCell'     'pDM' 'format'};
inpForm.defval = {range0  true       [0 0 0]     true       1     'plot'  };
inpForm.size   = {[-5 -6] [1 1]      [1 3]       [1 1]      [1 1] [1 -9]  };

inpForm.fname  = [inpForm.fname  {'angHeadAxis' 'lHeadAxis'     'rAxis'}];
inpForm.defval = [inpForm.defval {15             0.5           0.06    }];
inpForm.size   = [inpForm.size   {[1 1]          [1 1]         [1 1]   }];

inpForm.fname  = [inpForm.fname  {'pMultCell'    'pAxis'   'surfRes' 'rCoupling' 'sSpin'  'rSpin'}];
inpForm.defval = [inpForm.defval {false          true       30        0.05        1        0.06  }];
inpForm.size   = [inpForm.size   {[1 1]          [1 1]      [1 1]     [1 1]       [1 1]    [1 1] }];

inpForm.fname  = [inpForm.fname  {'angHeadSpin' 'lHeadSpin'     'aPlaneSpin'}];
inpForm.defval = [inpForm.defval {15             0.5           0.07         }];
inpForm.size   = [inpForm.size   {[1 1]          [1 1]         [1 1]        }];

inpForm.fname  = [inpForm.fname  {'labelAtom'  'cSpin' 'cCell' 'pCoupling' 'legend'}];
inpForm.defval = [inpForm.defval {true         'auto'    [0 0 0]    true    true    }];
inpForm.size   = [inpForm.size   {[1 1]        [1 -4]    [1 3]       [1 1]   [1 1]  }];

inpForm.fname  = [inpForm.fname  {'lineStyleCell' 'dAxis'       'dText' 'wSpace' 'rAtomData'}];
inpForm.defval = [inpForm.defval {'--'            [0.5;1.5;2.0] 0.2     10       false      }];
inpForm.size   = [inpForm.size   {[1 -1]          [3 1]         [1 1]   [1 1]    [1 1]      }];

inpForm.fname  = [inpForm.fname  {'rAtom' 'aEll' 'coplanar' 'fontSize' 'pNonMagAtom'}];
inpForm.defval = [inpForm.defval {0.3     0.3          0          12     true       }];
inpForm.size   = [inpForm.size   {[1 -8]  [1 1]        [1 1]      [1 1]  [1 1]      }];

inpForm.fname  = [inpForm.fname  {'rEll' 'eEll' 'pZeroCoupling' 'tooltip' 'cCoupling' 'cDM' 'cAtom'}];
inpForm.defval = [inpForm.defval {1          0.1          true             true      'auto' 'auto'   'auto' }];
inpForm.size   = [inpForm.size   {[1 1]      [1 1]        [1 1]            [1 1]     [1 -7] [1 -2]   [1 -3] }];

inpForm.fname  = [inpForm.fname  {'sEll' 'lwEll' 'dash' 'lineWidthCell' 'hFigure' 'hg'  'zoom' }];
inpForm.defval = [inpForm.defval {1      1       1      1                0        true  0      }];
inpForm.size   = [inpForm.size   {[1 1]  [1 1]   [1 1]  [1 1]            [1 1]    [1 1] [1 1]  }];

inpForm.fname  = [inpForm.fname  {'centerS' }];
inpForm.defval = [inpForm.defval {true      }];
inpForm.size   = [inpForm.size   {[1 1]     }];

param = sw_readparam(inpForm, varargin{:});

switch param.format
    case 'plot'
        plotmode = true;
    case 'jmol'
        plotmode = false;
        strOut = sprintf('BACKGROUND WHITE;\n');
    otherwise
        error('sw:plot:WrongInput','''format'' option has to be either ''plot'' or ''jmol''!');
end


lattice   = obj.lattice;
% The basis vectors in columns.
basisVector = obj.basisvector;

if plotmode
    % Handle of figure window to plot.
    if param.hFigure>0
        hFigure = param.hFigure;
        figure(hFigure);
    else
        hFigure = sw_getfighandle('sw_crystal');
    end
    
    if param.hg
        if isempty(hFigure)
            hFigure = sw_structfigure;
        end
    else
        hFigure = figure;
        if feature('usehg2')
            hFigure.Name = sprintf('Figure %d: SpinW : Crystal structure',hFigure.Number);
            hFigure.NumberTitle  = 'off';
            hFigure.DockControls = 'off';
            hFigure.PaperType    = 'A4';
            hFigure.Tag          = 'sw_crystal';
            %hFigure.Toolbar      = 'figure';
        else
            set(hFigure,...
                'Name',          sprintf('Figure %d: SpinW : Crystal structure',hFigure),...
                'NumberTitle',   'off',...
                'DockControls',  'off',...
                'PaperType',     'A4',...
                'Tag',           'sw_crystal',...
                'Toolbar',       'figure');
        end
        set(gcf,'color','w')
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
    end
    
    tooltip(hFigure,'<< Click on any object to get information! >>')
    
    cva = get(gca,'CameraViewAngle');
    [az, el] = view();
end

% change range, if the number of unit cells are given
if numel(param.range) == 3
    param.range = [ [0 0 0]' param.range(:)];
elseif numel(param.range) ~=6
    error('sw:plot:WrongInput','The plotting range is wrong, see doc sw.plot!');
end

%shift of the plot
wShift = -sum(basisVector * sum(param.range,2)/2,2);
wShift = repmat(wShift,[1 2]);
wShift = reshape(wShift',1,[]);

% Plot spins if spin variables exist and param.pSpin is true.
param.pSpin = param.pSpin && ~isempty(obj.mag_str.S);

%% save data to figure

if plotmode
    % delete previous object if a different one is going to be saved
    objOld = getappdata(hFigure,'obj');
    if isempty(objOld) || (objOld ~= obj)
        delete(objOld);
        setappdata(hFigure,'obj',copy(obj));
    end
    setappdata(hFigure,'param',param);
end

%% Calculated parameters

% basis vectors for unit cell
v0 = zeros(3,1);
v1 = basisVector(:,1);
v2 = basisVector(:,2);
v3 = basisVector(:,3);

%axis range to show
if plotmode
    param.wRange = [param.range(:,1)-param.wSpace param.range(:,2)+param.wSpace];
    axis(reshape((basisVector*param.wRange)',1,[])+wShift);
    zoom(10);
    
    %sphare surface for atoms
    [sp.x, sp.y, sp.z] = sphere(param.surfRes);
    % circles for the ellipsoid
    cr{1} = sw_circle([0 0 0]',[1 0 0]',1.01,param.surfRes);
    cr{2} = sw_circle([0 0 0]',[0 1 0]',1.01,param.surfRes);
    cr{3} = sw_circle([0 0 0]',[0 0 1]',1.01,param.surfRes);
    
    hold on
end

%% Plot the unit cell.


handle.unitCell = [];

if param.pCell
    idx = 1;
    
    path=[v0 v1 v1+v2 v2 v0 v3 v3+v1 v1 v3+v1 v3+v1+v2 v1+v2 v3+v1+v2 v3+v2 v2 v3+v2 v3];
    
    if param.pMultCell
        for ii = ceil(param.range(1,1)):floor(param.range(1,2)-1)
            for jj = ceil(param.range(2,1)):floor(param.range(2,2)-1)
                for kk = ceil(param.range(3,1)):floor(param.range(3,2)-1)
                    pathd = path+(v1*ii+v2*jj+v3*kk)*ones(1,16);
                    if plotmode
                        handle.unitCell(end+1) = plot3(pathd(1,:),pathd(2,:),pathd(3,:));
                    else
                        for ll = 2:size(pathd,2)
                            objid = sprintf('unitCell%d',idx);
                            idx = idx + 1;
                            strOut = [strOut jmol_command('line',objid,pathd(:,ll-1),pathd(:,ll),param.cCell)];
                        end
                    end
                end
            end
        end
    end
    if plotmode
        if isempty(handle.unitCell)
            handle.unitCell = plot3(path(1,:),path(2,:),path(3,:));
        end
    else
        if idx == 1
            for ll = 2:size(path,2)
                objid = sprintf('unitCell%d',idx);
                idx = idx + 1;
                strOut = [strOut jmol_command('line',objid,path(:,ll-1),path(:,ll),param.cCell)];
            end
        end
    end
    
    set(handle.unitCell,'Color',param.cCell,'LineStyle',param.lineStyleCell,'LineWidth',param.lineWidthCell);
    set(handle.unitCell,'Tag','unitCell');
    tooltip(handle.unitCell,['Unit cell \na=' sprintf('%5.3f',lattice.lat_const(1)) ' \nb=' sprintf('%5.3f',lattice.lat_const(2)) ' \nc=' sprintf('%5.3f',lattice.lat_const(3))])
end

%% Plot the coordinate system.
if param.pAxis
    
    r = -ones(3,1).*param.dAxis + v1*param.range(1,1) + v2*param.range(2,1) + v3*param.range(3,1);
    
    if plotmode
        handle.arrow = zeros(12,1);
        
        [handle.arrow(1:4)]  = sw_arrow(r,r+v1,param.rAxis,param.angHeadAxis,param.lHeadAxis,param.surfRes);
        [handle.arrow(5:8)]  = sw_arrow(r,r+v2,param.rAxis,param.angHeadAxis,param.lHeadAxis,param.surfRes);
        [handle.arrow(9:12)] = sw_arrow(r,r+v3,param.rAxis,param.angHeadAxis,param.lHeadAxis,param.surfRes);
        
        set(handle.arrow,'Facecolor',param.cAxis);
        set(handle.arrow,'Tag','arrow');
        
        v1n = v1/norm(v1)*param.dText;
        v2n = v2/norm(v2)*param.dText;
        v3n = v3/norm(v3)*param.dText;
        
        handle.arrowText(1) = text(r(1)+v1(1)+v1n(1),r(2)+v1(2)+v1n(2),r(3)+v1(3)+v1n(1),'a','fontSize',param.fontSize);
        handle.arrowText(2) = text(r(1)+v2(1)+v2n(1),r(2)+v2(2)+v2n(2),r(3)+v2(3)+v2n(2),'b','fontSize',param.fontSize);
        handle.arrowText(3) = text(r(1)+v3(1)+v3n(1),r(2)+v3(2)+v3n(2),r(3)+v3(3)+v3n(3),'c','fontSize',param.fontSize);
        
        set(handle.arrowText,'Tag','arrowText');
        set(handle.arrowText,'color',param.cAxis);
        tooltip(handle.arrow,'Arrow')
        tooltip(handle.arrowText,'Arrow text')
        
    else
        strOut = [strOut jmol_command('arrow','a-axis',param.rAxis,r,r+v1,param.cAxis)];
        strOut = [strOut jmol_command('arrow','b-axis',param.rAxis,r,r+v2,param.cAxis)];
        strOut = [strOut jmol_command('arrow','c-axis',param.rAxis,r,r+v3,param.cAxis)];
    end
end

if plotmode
    
    handle.atom        = [];
    handle.aniso       = [];
    handle.atomText    = [];
    handle.spinCircle  = [];
    handle.spinCircle2 = [];
    handle.spinArrow   = [];
    handle.coupling    = [];
    handle.DMcoupling  = [];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ATOMS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

atom       = obj.atom;
mAtom      = obj.matom;

nAtom      = size(atom.r,2);
nMagAtom   = size(mAtom.r,2);

atom.color = double(obj.unit_cell.color(:,atom.idx))/255;
atom.name  = obj.unit_cell.label(atom.idx);
atom.label = cell(1,nAtom);

% plot only the first word of every label
for ii = 1:nAtom
    labelTemp = strword(atom.name{ii},[1 2],true);
    atom.label{ii} = labelTemp{1};
    atom.name{ii}  = labelTemp{2};
    
    % labelTemp = strsplit(atom.name{ii});
    % wordIdx = find(cellfun(@numel,labelTemp));
    % atom.name{ii} = labelTemp{wordIdx(1)};
end

atom.rad   = ones(nAtom,1) * param.rAtom; % atomic radius

for ll = 1:nAtom
    % Saves atomic radius data if allowed and exists scaled with param.rAtom.
    if param.rAtomData && (sw_atomdata(atom.name{ll},'radius') > 0)
        if numel(param.rAtom)>1
            atom.rad(ll) = sw_atomdata(atom.name{ll},'radius')*param.rAtom(atom.idx(ll));
        else
            atom.rad(ll) = sw_atomdata(atom.name{ll},'radius')*param.rAtom;
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
% GENERATE HAMILTONIAN USING SYMMETRY
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[SS, SI] = intmatrix(obj,'plotmode',true,'extend',false,'zeroC',param.pZeroCoupling);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% COUPLINGS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

matrix   = obj.matrix;
coupling = obj.coupling;
nExt     = double(obj.mag_str.N_ext);

if isempty(SS.all)
    coupling.dl      = zeros(3,0);
    coupling.atom1   = zeros(1,0);
    coupling.atom2   = zeros(1,0);
    coupling.mat_idx = zeros(1,0);
    coupling.idx     = zeros(1,0);
    coupling.mat     = zeros(3,3,0);
    coupling.DM      = zeros(3,0);
else
    coupling.dl      = SS.all(1:3,:);
    coupling.atom1   = SS.all(4,:);
    coupling.atom2   = SS.all(5,:);
    coupling.mat_idx = SS.all(15,:);
    coupling.idx     = SS.all(16,:);
    coupling.mat     = reshape(SS.all(6:14,:),3,3,[]);
    coupling.DM      = (coupling.mat-permute(coupling.mat,[2 1 3]))/2;
    coupling.DM      = cat(2,coupling.DM(2,3,:),coupling.DM(3,1,:),coupling.DM(1,2,:));
end

nCoupling        = size(coupling.idx,2);

% dashes
dashList = zeros(1,numel(matrix.label));
for ii = 1:numel(matrix.label)
    dashList(ii) = matrix.label{ii}(end);
    if strcmp(char(dashList(ii)),'-')
        dashList(ii) = true;
    else
        dashList(ii) = false;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ANISOTROPY ELLIPSOIDS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

single_ion = obj.single_ion;
if numel(single_ion.aniso) ~= nMagAtom
    single_ion.aniso = zeros(1,nMagAtom);
end
if numel(single_ion.g) ~= nMagAtom
    single_ion.g = zeros(1,nMagAtom);
end

% select anisotropy matrices or g-tensor
switch param.sEll
    case 1
        single_ion.mat = SI.aniso;
        single_ion.idx = single_ion.aniso;
    case 2
        single_ion.mat = SI.g;
        single_ion.idx = single_ion.g;
    case 0
    otherwise
        error('sw:plot:WrongParameter','The sEll parameter has to be 0,1 or 2!');
end

single_ion.ell   = zeros(3,3,nMagAtom);
single_ion.label = cell(1,nMagAtom);

if param.sEll>0
    % Creates the anisotropy/g-tensor ellipsoids.
    for ll = 1:nMagAtom
        % Selects the matrix index to plot
        cIdx = single_ion.idx(ll);
        
        if any(cIdx)
            
            % Calculates the main radiuses of the ellipsoid.
            %[V, R] = eig(single_ion.mat(:,:,ll));
            [V, R0] = eigorth(single_ion.mat(:,:,ll),1e-5);
            % Creates positive definite matrix by adding constant to all
            % eigenvalues.
            %R0      = diag(R);
            epsilon = sqrt(sum(R0.^2))*param.eEll;
            dR      = 1./(R0-min(R0)+epsilon);
            dR0     = max(dR);
            if dR0 ~= 0
                dR = dR/dR0;
            else
                dR = [1 1 1];
            end
            
            R = diag(dR)*param.rEll;
            
            single_ion.ell(:,:,ll) = V*R;
            single_ion.label{ll}   = matrix.label{cIdx};
        end
    end
end

if plotmode
    
else
    idxe = 1;
    idxs = 1;
    idxa = 1;
    idxc = 1;
end

%% Plot all atoms and couplings in selected range.
for ii = floor(param.range(1,1)):floor(param.range(1,2))
    for jj = floor(param.range(2,1)):floor(param.range(2,2))
        for kk = floor(param.range(3,1)):floor(param.range(3,2))
            dCell = [ii;jj;kk];
            % Store the index of magAtom
            llMagAtom = 1;
            for ll = 1:nAtom
                % color of the atoms
                if strcmpi(param.cAtom,'auto')
                    AColor = atom.color(:,ll);
                else
                    AColor = param.cAtom/255;
                end
                if strcmpi(param.cSpin,'auto')
                    MColor = atom.color(:,ll);
                else
                    MColor = param.cSpin/255;
                end
                
                % Index of spin in the extended unit cell
                llSpin = mod(dCell,nExt')'*[1; nExt(1); nExt(1)*nExt(2)]*nMagAtom;
                
                % Atomic position in lattice units.
                rLat = atom.r(:,ll) + dCell;
                
                % Plot if the atom is within the plotting range.
                if all((rLat <= param.range(:,2)) & (rLat >= param.range(:,1)))
                    rPlot = basisVector*rLat;
                    
                    % Draw anisotropy semi-transparent ellipsoid.
                    if atom.mag(ll) && param.sEll>0
                        aIdx = single_ion.idx(llMagAtom);
                        if any(aIdx) && (norm(obj.matrix.mat(:,:,aIdx))>0)
                            
                            if plotmode
                                ell.xyz = single_ion.ell(:,:,llMagAtom)*[sp.x(:) sp.y(:) sp.z(:)]';
                                
                                ell.x = reshape(ell.xyz(1,:),[1 1]*param.surfRes+1);
                                ell.y = reshape(ell.xyz(2,:),[1 1]*param.surfRes+1);
                                ell.z = reshape(ell.xyz(3,:),[1 1]*param.surfRes+1);
                                
                                sAniso = surf(ell.x+rPlot(1),ell.y+rPlot(2),ell.z+rPlot(3));
                                
                                handle.aniso(atom.idx(ll),end+1) = sAniso;
                                set(sAniso,'LineStyle','none');
                                set(sAniso,'FaceAlpha',param.aEll);
                                set(sAniso,'FaceColor',double(obj.matrix.color(:,aIdx))/255);
                                set(sAniso,'Tag',['aniso_' atom.label{ll}]);
                                
                                tooltip(sAniso,[single_ion.label{llMagAtom} ' anisotropy\n' strmat(single_ion.mat(:,:,llMagAtom))]);
                                
                                if param.lwEll > 0
                                    % draw the main circles of the ellipsoids
                                    ell.c1 = single_ion.ell(:,:,llMagAtom)*cr{1};
                                    ell.c2 = single_ion.ell(:,:,llMagAtom)*cr{2};
                                    ell.c3 = single_ion.ell(:,:,llMagAtom)*cr{3};
                                    
                                    sC(1) = plot3(ell.c1(1,:)+rPlot(1),ell.c1(2,:)+rPlot(2),ell.c1(3,:)+rPlot(3));
                                    sC(2) = plot3(ell.c2(1,:)+rPlot(1),ell.c2(2,:)+rPlot(2),ell.c2(3,:)+rPlot(3));
                                    sC(3) = plot3(ell.c3(1,:)+rPlot(1),ell.c3(2,:)+rPlot(2),ell.c3(3,:)+rPlot(3));
                                    set(sC,'LineWidth',param.lwEll);
                                    set(sC,'Color',[0 0 0]);
                                    set(sC,'Tag',['aniso_circle_' atom.label{ll}]);
                                    handle.aniso(atom.idx(ll),end+(1:3)) = sC;
                                end
                                
                            else
                                objid = sprintf('anisotropy%d',idxe);
                                idxe = idxe + 1;
                                jmol_command('ellipsoid',objid,rPlot,ell.xyz(:,1),ell.xyz(:,2),ell.xyz(:,3),obj.matrix.color(:,aIdx));
                            end
                            
                        end
                    end
                    
                    if (atom.mag(ll) || param.pNonMagAtom) && (param.rAtom(1)>0)
                        % Plot the label of the atom.
                        if plotmode
                            if (~any(dCell)) && param.labelAtom
                                
                                handle.atomText(atom.idx(ll),end+1) = text(...
                                    'Position',             rPlot+param.dText+atom.rad(ll)/sqrt(3),...
                                    'String',               sprintf('%s(%d)_{%d}',atom.label{ll},atom.idx(ll),ll),...
                                    'VerticalAlignment',    'bottom',...
                                    'Tag',                  ['atomText_' atom.label{ll} ' '],...
                                    'fontSize',             param.fontSize);
                                tooltip(handle.atomText(atom.idx(ll),end),[atom.label{ll} ' atom \nUnit cell: \n' sprintf('[%d, %d, %d]',dCell) '\nAtomic position: \n' sprintf('[%6.3f, %6.3f, %6.3f] ',atom.r(:,ll))]);
                            end
                            
                            % Creates the sphere of the atom.
                            aSphere = surf(sp.x*atom.rad(ll)+rPlot(1), sp.y*atom.rad(ll)+rPlot(2), sp.z*atom.rad(ll)+rPlot(3));
                            set(aSphere,'Tag',['atom_' atom.label{ll}]);
                            handle.atom(atom.idx(ll),end+1) = aSphere;
                            set(aSphere,'LineStyle','none');
                            set(aSphere,'FaceColor',AColor);
                            
                            tooltip(aSphere,[atom.name{ll} ' atom (' atom.label{ll} ') \nUnit cell: \n' sprintf('[%d, %d, %d]',dCell) '\nAtomic position: \n' sprintf('[%6.3f, %6.3f, %6.3f] ',atom.r(:,ll))]);
                        else
                            objid = sprintf('%s(%d)_{%d}%d',atom.label{ll},atom.idx(ll),ll,idxa);
                            idxa = idxa + 1;
                            strOut = [strOut jmol_command('sphere',objid,atom.rad(ll),rPlot,AColor*255)]; %#ok<*AGROW>
                        end
                    end
                    
                    % Plots magnetic moments.
                    if param.pSpin && atom.mag(ll)
                        % translation in lat.units to calculate moment
                        transl  = floor(dCell./nExt').*nExt';
                        % angles of rotation of the moment
                        phi   = obj.mag_str.k*transl*2*pi;
                        selS  = obj.mag_str.S(:,llMagAtom+llSpin);
                        plotS = param.sSpin*sw_rot(obj.mag_str.n,phi,selS);
                        %plotS   = param.sSpin*obj.mag_str.S(:,llMagAtom+llSpin);
                        lengthS = norm(obj.mag_str.S(:,llMagAtom+llSpin));
                        
                        if obj.symbolic
                            symVar1 = symvar(plotS);
                            symVar2 = symvar(lengthS);
                            if ~isempty(symVar1)
                                plotS = double(subs(plotS,symVar1,ones(1,numel(symVar1))));
                            end
                            if ~isempty(symVar2)
                                lengthS = double(subs(lengthS,symVar2,ones(1,numel(symVar2))));
                            end

                            %plotS = double(subs(plotS,1));
                            %lengthS = double(subs(lengthS,1));
                        end
                        
                        if param.coplanar
                            % TODO
                            % Circle for planar structure and maybe cones :)
                            if plotmode
                                cPoints = sw_circle(rPlot,nSpin,norm(plotS),param.surfRes);
                                C1 = plot3(cPoints(1,:),cPoints(2,:),cPoints(3,:));
                                handle.spinCircle(atom.idx(ll),end+1) = C1;
                                set(C1,'Color',MColor);
                                set(C1,'Tag','spinPlane');
                                tooltip(C1,['Plane of the ordered moment: \n' sprintf('[%6.3f, %6.3f, %6.3f]',nSpin)]);
                                
                                C2 = sw_circlesurf(rPlot,nSpin,norm(plotS),param.surfRes);
                                handle.spinCircle2(atom.idx(ll),end+1) = C2;
                                set(C2,'LineStyle','none');
                                set(C2,'FaceAlpha',param.aPlaneSpin);
                                set(C2,'FaceColor',MColor);
                                tooltip(C2,['Plane of the ordered moment: \n' sprintf('[%6.3f, %6.3f, %6.3f]',nSpin)]);
                            end
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
                        
                        if plotmode
                            if param.centerS
                                hArrow  = sw_arrow(rPlot-plotS,rPlot+plotS,param.rSpin,param.angHeadSpin,param.lHeadSpin,param.surfRes);
                            else
                                hArrow  = sw_arrow(rPlot,rPlot+plotS,param.rSpin,param.angHeadSpin,param.lHeadSpin,param.surfRes);
                            end
                            set(hArrow,'Tag','spinArrow');
                            % save information into the spin arrow
                            % [a b c idx] unit cell position and spin index
                            setappdata(hArrow(1),'dat',[dCell' llMagAtom+llSpin]);
                            handle.spinArrow(atom.idx(ll),end+(1:numel(hArrow))) = hArrow;
                            set(hArrow,'FaceColor',MColor);
                            set(hArrow,'LineStyle','none');
                            if intS
                                if halfS
                                    str0 = ['S=' sprintf('%d',lengthS2) '/2'];
                                    %tooltip(hArrow,['S=' sprintf('%d',lengthS2) '/2 magnetic moment \nxyz components:  \n' sprintf('[%6.3f, %6.3f, %6.3f]',double(subs(obj.mag_str.S(:,llMagAtom+llSpin),1)))])
                                else
                                    str0 = ['S=' sprintf('%d',lengthS)];
                                    %tooltip(hArrow,['S=' sprintf('%d',lengthS) ' magnetic moment \nxyz components:  \n' sprintf('[%6.3f, %6.3f, %6.3f]',double(subs(obj.mag_str.S(:,llMagAtom+llSpin),1)))])
                                end
                            else
                                str0 = ['S=' sprintf('%4.2f',lengthS)];
                                %tooltip(hArrow,['S=' sprintf('%4.2f',lengthS) ' magnetic moment \nxyz components:  \n' sprintf('[%6.3f, %6.3f, %6.3f]',double(subs(obj.mag_str.S(:,llMagAtom+llSpin),1)))])
                            end
                            if obj.symbolic
                                symVar1 = symvar(obj.mag_str.S(:,llMagAtom+llSpin));
                                
                                if ~isempty(symVar1)
                                    spin0 = double(subs(obj.mag_str.S(:,llMagAtom+llSpin),symVar1,ones(1,numel(symVar1))));
                                end

                                str1 = sprintf('[%6.3f, %6.3f, %6.3f]',double(spin0));
                            else
                                str1 = sprintf('[%6.3f, %6.3f, %6.3f]',double(obj.mag_str.S(:,llMagAtom+llSpin)));
                            end
                            tooltip(hArrow,[str0 ' magnetic moment \nxyz components:  \n' str1]);
                        else
                            objid = sprintf('spinArrow%d',idxs);
                            idxs = idxs + 1;
                            if param.centerS
                                strOut = [strOut jmol_command('arrow',objid,param.rSpin,rPlot-plotS,rPlot+plotS,MColor*255)];
                            else
                                strOut = [strOut jmol_command('arrow',objid,param.rSpin,rPlot,rPlot+plotS,MColor*255)];
                            end
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
                if coupling.mat_idx(ll)
                    rLat1 = mAtom.r(:,coupling.atom1(ll)) + dCell;
                    rLat2 = mAtom.r(:,coupling.atom2(ll)) + double(coupling.dl(:,ll)) + dCell;
                    
                    if all((rLat1<=param.range(:,2))&(rLat1>=param.range(:,1))&(rLat2<=param.range(:,2))&(rLat2>=param.range(:,1)))
                        rPlot1    = basisVector*rLat1;
                        rPlot2    = basisVector*rLat2;
                        if strcmpi(param.cCoupling,'auto')
                            cColor    = double(matrix.color(:,coupling.mat_idx(ll)))/255;
                        else
                            cColor    = param.cCoupling/255;
                        end
                        if strcmpi(param.cDM,'auto')
                            DMColor    = double(matrix.color(:,coupling.mat_idx(ll)))/255;
                        else
                            DMColor    = param.cDM/255;
                        end
                        
                        % plot normal cylinder
                        if param.pCoupling
                            if plotmode
                                dashS = (param.dash*dashList(coupling.mat_idx(ll)));
                                if dashS>0
                                    cCylinder = sw_cylinder(rPlot1,rPlot2,param.rCoupling*1.2,param.surfRes,dashS);
                                else
                                    cCylinder = sw_cylinder(rPlot1,rPlot2,param.rCoupling,param.surfRes,0);
                                end
                                handle.coupling(coupling.idx(ll),end+(1:numel(cCylinder))) = cCylinder;
                                set(cCylinder,'FaceColor',cColor);
                                set(cCylinder,'Tag',['coupling_' matrix.label{coupling.mat_idx(ll)}]);
                                
                                tooltip(cCylinder,[matrix.label{coupling.mat_idx(ll)} ' coupling\nValue:\n' strmat(matrix.mat(:,:,coupling.mat_idx(ll)))]);
                            else
                                objid = sprintf('coupling_%s%d',matrix.label{coupling.mat_idx(ll)},idxc);
                                idxc = idxc + 1;
                                strOut = [strOut jmol_command('cylinder',objid,param.rCoupling,rPlot1,rPlot2,cColor*255)];
                            end
                        end
                        % plot DM vector
                        if param.pDM>0 && (norm(coupling.DM(1,:,ll))>1e-5)
                            rCent = (rPlot1+rPlot2)/2;
                            vDM   = coupling.DM(1,:,ll)'*param.pDM;
                            if plotmode
                                hDM   = sw_arrow(rCent,rCent+vDM,param.rSpin,param.angHeadSpin,param.lHeadSpin,param.surfRes);
                                handle.DMcoupling(coupling.idx(ll),end+(1:4)) = hDM;
                                set(hDM,'FaceColor',DMColor);
                                set(hDM,'Tag',['coupling_' matrix.label{coupling.mat_idx(ll)}]);
                                
                                tooltip(hDM,[matrix.label{coupling.mat_idx(ll)} ' DM coupling\nValue:\n' strmat(vDM')]);
                            else
                                objid = sprintf('coupling_%s%d',matrix.label{coupling.mat_idx(ll)},idxc);
                                idxc = idxc + 1;
                                strOut = [strOut jmol_command('arrow',objid,param.rSpin,rCent,rCent+vDM,DMColor*255)];
                                
                            end
                            if ~param.pCoupling
                                
                                if plotmode
                                    rPlot = [rPlot1 rPlot2];
                                    hLine = line(rPlot(1,:),rPlot(2,:),rPlot(3,:),'LineStyle','--','LineWidth',2,'Color',[0 0 0]);
                                    set(hLine,'Tag',['coupling_' matrix.label{coupling.mat_idx(ll)}]);
                                    handle.DMcoupling(coupling.idx(ll),end+1) = hLine;
                                else
                                    objid = sprintf('coupling_%s%d',matrix.label{coupling.mat_idx(ll)},idxc);
                                    idxc = idxc + 1;
                                    strOut = [strOut jmol_command('line',objid,rPlot1,rPlot2,[0 0 0])];
                                end
                            end
                        end
                    end
                end
                
            end
        end
    end
end

%% for Jmol output quit here
if ~plotmode
    strOut = [strOut sprintf('set drawHover on;')];
    varargout{1} = strOut;
    return;
end

%% Do the rest.

view([az el]);
hold off

if isappdata(hFigure,'handle')
    oldHandle = getappdata(hFigure,'handle');
    hName = fieldnames(oldHandle);
    
    for ii = 1:length(hName)
        h0 = reshape(oldHandle.(hName{ii}),1,[]);
        h0(h0 == 0) = [];
        h0(~ishandle(h0)) = [];
        delete(h0);
    end
end

% put all objects into a hgtransform object for rotation
if param.hg
    h  = getappdata(hFigure,'h');
    h2 = hgtransform('Parent',h);
    hName = fieldnames(handle);
    
    for ii = 1:length(hName)
        h0 = reshape(handle.(hName{ii}),1,[]);
        h0(h0 == 0) = [];
        h0(~ishandle(h0)) = [];
        set(h0,'Parent',h2);
        set(h0,'Clipping','Off');
    end
    T = makehgtform('translate',-sum(basisVector * sum(param.range,2)/2,2)');
    set(h2,'Matrix',T);
    
end

set(gca,'CameraViewAngle',cva);
handle.light = camlight('right');
set(handle.light,'Tag','light');

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
    handle.lRect = [];
    for ii = 1:length(matrix.label)
        handle.lRect(end+1) = rectangle('Position',[5 lHeight-20*ii+6 sRadius*2 sRadius]);
        if dashList(ii)
            % legend for dashed matrices
            handle.lRect(end+1) = rectangle('Position',[6 lHeight-20*ii+6.9 sRadius*2/3-1 sRadius-1.5]);
            set(handle.lRect(end),'FaceColor',double(matrix.color(:,ii))/255);
            set(handle.lRect(end),'EdgeColor',double(matrix.color(:,ii))/255);
            
            handle.lRect(end+1) = rectangle('Position',[5+sRadius*2/3*2 lHeight-20*ii+6.9 sRadius*2/3-1 sRadius-1.5]);
            set(handle.lRect(end),'FaceColor',double(matrix.color(:,ii))/255);
            set(handle.lRect(end),'EdgeColor',double(matrix.color(:,ii))/255);
        else
            set(handle.lRect(end),'FaceColor',double(matrix.color(:,ii))/255);
        end
        
        if dashList(ii)
            handle.lText(ii) = text(30,(lHeight-20*ii+10),matrix.label{ii}(1:end-1),...
                'fontSize',param.fontSize);
        else
            handle.lText(ii) = text(30,(lHeight-20*ii+10),matrix.label{ii},...
                'fontSize',param.fontSize);
        end
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

axes(cAxis);

%% Stuff...

setappdata(hFigure,'handle',handle);

if get(gca,'CameraViewAngle') == 0.6
    camva('auto');
end;

% apply zoom setting
cva = get(gca,'CameraViewAngle');
set(gca,'CameraViewAngle',cva*(1.5)^(-param.zoom));


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
    string = [string '|'];
    for jj = 1:size(mat,2)
        if isa(mat,'sym')
            string = [string char(mat(ii,jj))];
        else
            string = [string sprintf(['%' num2str(p1) '.' num2str(p2) 'f '],mat(ii,jj))];
        end
    end
    string = [string '|\n'];
end

end

function tooltip(hObject, string)
% Creates tooltip for graphical objects, with string.

set(hObject,'buttondownfcn',['oTemp=findobj(''tag'',''tooltip''); if ~isempty(oTemp) set(oTemp,''string'',sprintf(''' string '''));end; clear oTemp;']);

end

function str = jmol_command(shape, varargin)
% create Jmol script lines
% shape:
%   cylinder    (id,radius,R1,R2,color)
%   sphere      (id,radius,R1,color)
%   ellipsoid   (id,R1,ax1,ax2,ax3,color)
%   arrow       (id,radius,R1,R2,color)
%   line        (id,R1,R2,color)
%

switch shape
    case 'cylinder'
        str = sprintf('DRAW ID "%s" CYLINDER RADIUS %5.3f {%5.3f,%5.3f,%5.3f} {%5.3f,%5.3f,%5.3f} COLOR {%3d,%3d,%3d} "%s";\n',varargin{1},varargin{2},varargin{3},varargin{4},varargin{5},varargin{1});
    case 'sphere'
        R = varargin{2};
        str = sprintf('ELLIPSOID ID "%s" CENTER {%5.3f,%5.3f,%5.3f} AXES {%5.3f,0,0} {0,%5.3f,0} {0,0,%5.3f} COLOR {%3d,%3d,%3d};\n',varargin{1},varargin{3},R,R,R,varargin{4});
    case 'ellipsoid'
        str = sprintf('ELLIPSOID ID "%s" CENTER {%5.3f,%5.3f,%5.3f} AXES {%5.3f,%5.3f,%5.3f} {%5.3f,%5.3f,%5.3f} {%5.3f,%5.3f,%5.3f} COLOR {%3d,%3d,%3d};\n',...
            varargin{1},varargin{2},varargin{3},varargin{4},varargin{5},varargin{6});
    case 'arrow'
        str = sprintf('DRAW ID "%s" ARROW RADIUS %5.3f {%5.3f,%5.3f,%5.3f} {%5.3f,%5.3f,%5.3f} COLOR {%3d,%3d,%3d} "%s";\n',varargin{1},varargin{2},varargin{3},varargin{4},varargin{5},varargin{1});
    case 'line'
        str = sprintf('DRAW ID "%s" LINE RADIUS 0.01 {%5.3f,%5.3f,%5.3f} {%5.3f,%5.3f,%5.3f} COLOR {%3d,%3d,%3d} "%s";\n',varargin{1},varargin{2},varargin{3},varargin{4},varargin{1});
    otherwise
        error('sw:plot:Internal','Internal error in sw.plot(), please submit a bug report!');
end


end