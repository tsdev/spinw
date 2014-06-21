function sw_animate(spectra, varargin)
% animates normal magnon modes
%
% SW_ANIMATE(spectra, 'option1', value1 ...)
%
% Animates the selected normal magnon mode on a structure plot. The
% animation will run until the given number of full rotations or until it
% is aborted using the Ctrl+C command in the Matlab Command Window.
%
% Options:
%
% pAmp      Amplitude of the motion, default is 1.
% rps       Rotation per second, default is 1.
% fps       Frame per second, default is 25 for continuous rotation.
% phi0      Starting phase in radian, default is zero.
% nCycle    Number of full rotations, use Inf for infinite number
%           of rotations (stops with Ctrl+C). Default is Inf. For zero only
%           the t=0 phase is drawn.
% tMax      Maximum time to run in seconds. Default is zero when nCycle
%           determines the running time.
% dMagnon   Thickness of the spin precession ellipse, if zero it won't be
%           plotted. Default is 0.05.
% cMagnon   Color of the ellipses of the spin precession, [R G B] color
%           code, where R, G, B = 0..255. Defeult is 'auto', when the color
%           of the ellipses are the same as the magnetic atoms.
% aMagnon   Alpha value of the ellipses, default is 0.5 for semi
%           transparent ellipses.
% surfRes   Number of points on the surface mesh, default is 30.
%
% Example:
%
% spectra = crystal.spinwave({[0 0 0] [1 0 0] 20},'saveT',true);
% spectra = sw_spinmotion(spectra,'iMagnon',1,'Q',[1/2 0 0])
% plot(crystal)
% sw_animate(spectra)
%
% The above example will animate the normal magnon mode with index at
% Q = (1/2 0 0).
%
% See also SW.SPINWAVE, SW_SPINMOTION.
%

if nargin == 0
    help sw_animate;
    return;
end

% sw object of the model
obj = spectra.obj;

%persistent isrunning

inpForm.fname  = {'hFigure' 'pAmp' 'rps' 'fps' 'nCycle' 'phi0' 'tMax' 'dMagnon' 'cMagnon' 'aMagnon' 'surfRes'};
inpForm.defval = {0         1      1     25    Inf      0      0      0.05      'auto'    0.5       30       };
inpForm.size   = {[1 1]     [1 1]  [1 1] [1 1] [1 1]   [1 1]   [1 1]  [1 1]     [1 -1]    [1 1]     [1 1]    };

param = sw_readparam(inpForm, varargin{:});

% figure handle
if param.hFigure>0
    hFigure = param.hFigure;
    figure(hFigure);
else
    hFigure = sw_getfighandle('sw_crystal');
end

if isempty(hFigure)
    error('sw_animate:WrongInput','Couln''t find plot crystal structure plot window, use sw.plot() function!');
end

% struct storing all object handles in plot window
pHandle = getappdata(hFigure,'handle');

% parent hgtransform of all objects in plot window
h = get(pHandle.spinArrow(1),'Parent');

% number of spin arrows on plot
nVector = numel(pHandle.spinArrow)/4;

if nVector == 0
    warning('sw_animate:NothingToPlot','No spin arrows on plot to animate!');
    return;
end

% store handles of the hgtransform objects
hgMagnon = zeros(1,nVector);

% store data of the spin vectors
vectDat = zeros(nVector,4);

% loop through every spin arrow and insert a hgtransform
for ii = 1:nVector
    hgMagnon(ii) = hgtransform('Parent',h);
    set(pHandle.spinArrow((ii-1)*4+1:ii*4),'Parent',hgMagnon(ii));
    
    % we need the position and index of the magnetic moment stored in the
    % plotted arrows
    vectDat(ii,:) = getappdata(pHandle.spinArrow(ii*4-3),'dat');
end

% define the cleanup routine: plot original moment directions
cleanupObj = onCleanup(@()cleanMagnon(hgMagnon, pHandle.spinArrow));

% position of the spin for translation
posVect = obj.basisvector*vectDat(:,1:3)';

% classical vector directions of the magnetic moments
M0 = obj.magtable;

% for static plot only 1 phi angle is plotted
if param.nCycle == 0
    nPhi = 1;
else
    nPhi = param.fps/param.rps;
end

% vector of rotation angles around the equilibrium with phi0 starting angle
% at t=0 time.
vPhi = linspace(0,2*pi,nPhi+1);
vPhi = vPhi(1:end-1)+param.phi0;

iCycle = 0;

nCycle = param.nCycle;
if nCycle == 0
    nCycle = 1;
end

if ~isfield(spectra,'motion')
    error('sw_animate:WrongInput','Use sw_spinmotion() function first to select magnon mode to plot!');
end

% phase difference from the magnon momentum
phiK = spectra.hkl(:,spectra.motion.iQ)'*vectDat(:,1:3)'*2*pi;

if param.tMax > 0
    tic;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DRAW ELLIPSE OF MOTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% delete previous ellipses
if isfield(pHandle,'magnonEll')
    delete(pHandle.magnonEll);
end
pHandle.magnonEll = [];

if param.dMagnon > 0
    % major and minor axes of the spin precession ellipses
    Amp1 = spectra.motion.Amp1;
    Amp2 = spectra.motion.Amp2;
    
    cMagnon1 = M0.M;
    cMagnon2 = M0.M + bsxfun(@mrdivide,M0.M,sqrt(sum(M0.M.^2,1)))*param.dMagnon;
    
    % draw ellipse
    % TODO rotation for incommensurate magnetic structures
    for ii = 1:nVector
        idx = vectDat(ii,4);
        cEll1 = cMagnon1(:,idx) + posVect(:,ii);
        cEll2 = cMagnon2(:,idx) + posVect(:,ii);
        pHandle.magnonEll(end+(1:3)) = sw_ellcylinder(cEll1,cEll2,[1 1],param.surfRes,param.pAmp*Amp1(:,idx),param.pAmp*Amp2(:,idx), true);
    end
    set(pHandle.magnonEll,'tag','magnonEll');
    
    % set ellipse color
    if strcmp(param.cMagnon,'auto')
        for ii = 1:nVector
            cEll = obj.unit_cell.color(:,M0.atom(vectDat(ii,4)));
            set(pHandle.magnonEll,'FaceColor',double(cEll)/255);
        end
    else
        set(pHandle.magnonEll,'FaceColor',double(param.cMagnon)/255);
    end
    set(pHandle.magnonEll,'Clipping','Off');
    set(pHandle.magnonEll,'FaceAlpha',param.aMagnon);
    set(pHandle.magnonEll,'parent',h);
    setappdata(hFigure,'handle',pHandle);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ANIMATE!
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
while iCycle < nCycle
    for tt = 1:nPhi
        phit = vPhi(tt);
        for ii = 1:nVector
            % add a phase due to the magnon momentum
            phi = phit + phiK(ii);
            
            idx = vectDat(ii,4);
            selS = M0.M(:,idx);
            
            % translate vector to the origin
            M1 = eye(4);
            M1(1:3,4) = -posVect(:,ii);
            
            % move the spin out of the equilibrium along Amp1
            Avect1 = Amp1(:,idx);
            Avect2 = Amp2(:,idx);
            
            % minor major ellipse axes
            Amp = [norm(Avect1) norm(Avect2)];
            
            % calculated amplitude for the magnon
            Acalc = prod(Amp)/sqrt((Amp(2)*cos(phi))^2+(Amp(1)*sin(phi))^2);
            
            % move spin out of equilibrium
            %M2 = makehgtform('axisrotate',Avect2,Acalc/norm(selS)*2*pi*param.pAmp);
            M2 = makehgtform('axisrotate',Avect2,asin(Acalc*param.pAmp/norm(selS)));
            
            % rotate around the spin direction
            M3 = makehgtform('axisrotate',selS,phi);
            
            set(hgMagnon(ii),'Matrix',inv(M1)*M3*M2*M1); %#ok<MINV>
        end
        pause(1/param.fps);
        if (param.tMax > 0) && (toc>param.tMax)
            return;
        end
    end
    
    iCycle = iCycle + 1;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% STATIC SPIN DEVIATION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if param.nCycle == 0
    while 1
        if (param.tMax > 0) && (toc>param.tMax)
            return
        end
        pause(1/param.fps);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PUT SPINS BACK TO ORIGINAL DIRECTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function cleanMagnon(hgMagnon, hSpin)
        % remove magnon hgtransform objects
        
        hParent = get(hgMagnon(1),'Parent');
        set(hSpin,'Parent',hParent);
        delete(hgMagnon);
    end

end