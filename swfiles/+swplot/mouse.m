function mouse(hFigure, perspective)
% adds mouse callbacks to swplot figure
%
% SWPLOT.MOUSE({hFigure},{perspective})
%
% Adds rotation and zoom functionality to any swplot figure. The following
% mouse actions are supported:
%   - mouse-drag        Rotation of objects.
%   - ctrl+mouse-drag   Shift of objects (pan).
%   - mouse-wheel       Zoom of objects.
%   - ctrl+mouse-wheel  Change perspective and switch to perspective
%                       projection.
%
% Input:
%
% hFigure       Handle of the swplot figure. Default is the selected
%               figure.
% perspective   String determines whether camera projection mode is changed
%               automatically between orthographic (zooming withouth CTRL 
%               key pressed) and perspective (zooming with CTRL key
%               pressed):
%                   'auto'      Automatic switching (default).
%                   'fix'       No switching.
%
% See also CAMPROJ.
%

% Using code from:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright (c) 2008 Andrea Tagliasacchi
% All Rights Reserved
% email: ata2 at cs dot nospam dot sfu dot ca
% $Revision: 1.0$ 10 February 2008
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin == 0
    hFigure = swplot.activefigure;
end

if nargin < 2
    perspective = 'auto';
end

% automatic perspective/orthographic switch on zooming
switch perspective
    case 'auto'
        AP = true;
    case 'fix'
        AP = false;
    otherwise
        error('mouse:WrongInput','Value of perspective is invalid!');
end

hAxis = getappdata(hFigure,'axis');

% get hgtransform
hTranslate = getappdata(hFigure,'h');

if ~isempty(hTranslate)
    % assign mouse rotation
    hRotate    = get(hTranslate,'Parent');
    
    % create hgtransform only if requested
    % take care of fine object rotations via the mouse
    set(hFigure,'WindowButtonMotionFcn',@motion_callback);
    set(hFigure,'WindowButtonDownFcn',  @buttondown_callback);
    set(hFigure,'WindowButtonUpFcn',    @buttonup_callback);
    
    % create global variables for callback communication
    mousestatus = 'buttonup';
    START       = [0 0 0];
    MR0          = get(hRotate,'Matrix');
    MT0          = get(hTranslate,'Matrix');
    CTRL0       = false;
end

% zoom function
set(hFigure,'WindowScrollWheelFcn', @wheel_callback);

    function wheel_callback(~, event)
        % zoom in/out for event
        cva = get(hAxis,'CameraViewAngle');
        scale = 1.02^(event.VerticalScrollCount);
        % is ctrl pressed
        ctrl = ismember('control',get(gcbf,'currentmodifier'));
        
        if ctrl
            % change to perspective view
            if AP && strcmp(camproj,'orthographic')
                camproj('perspective');
            end
            % perspective change, but don't get too close
            newCva = sind(cva/2)/scale;
            %cPos   = get(hAxis,'CameraPosition');
            %fprintf('%g %g\n',[newCva cPos(3)]);
            
            %if cPos(3)*scale>1e4
            %    set(hAxis,'CameraPosition',[0 0 1e4]);
            %elseif newCva < 0.5
            set(hAxis,'CameraPosition',get(hAxis,'CameraPosition')*scale);
            set(hAxis,'CameraViewAngle',2*asind(newCva));
            %end
        else
            % change to orthographic view
            if AP && strcmp(camproj,'perspective')
                camproj('orthographic');
            end

            % standard Matlab call when java fails
            set(hAxis,'CameraViewAngle',cva*scale);
        end
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright (c) 2008 Andrea Tagliasacchi
% All Rights Reserved
% email: ata2 at cs dot nospam dot sfu dot ca
% $Revision: 1.0$ 10 February 2008
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% motion callback, event "every" mouse movement
    function motion_callback(src, ~)
        
        if strcmp(mousestatus, 'buttonup')
            % do nothing
            return
        end
        
        % is ctrl key pressed?
        ctrl = ismember('control',get(gcbf,'currentmodifier'));
        
        if ctrl~=CTRL0
            % reset
            buttondown_callback(src, []);
            return
        end
        
        % retrieve the current point
        P = get(gcf,'CurrentPoint');
        % retrieve window geometry
        HGEOM = get(src, 'Position');
        width  = HGEOM(3);
        height = HGEOM(4);

        if ctrl
            % pan
            % camera distance
            dCam = norm(get(hAxis,'CameraPosition'));
            % camera view angle
            cva = get(hAxis,'CameraViewAngle');
            % scale from window plane to object plane
            scale = dCam*sind(cva/2)/min([width height]/2);
            dP = [scale.*P 0] - START;
            MT = MT0;
            MT(1:3,4) = MT(1:3,4)+inv(MR0(1:3,1:3))*dP(:); %#ok<MINV>
            set(hTranslate,'Matrix',MT);
        else
            % rotate
            % transform in sphere coordinates (3,4) = (WIDTH, HEIGHT)
            P = point_on_sphere( P, width, height );
            %%%%% ARCBALL COMPUTATION %%%%%
            %if strcmp(mousestatus, 'buttondown')
            % compute angle and rotation axis
            rot_dir = cross( START, P ); rot_dir = rot_dir / norm( rot_dir );
            rot_ang = acos( dot( P, START ) );
            
            % convert direction in model coordinate system
            rot_dir = MR0\[rot_dir,0]';
            rot_dir = rot_dir(1:3);
            if norm(rot_dir) > 0
                rot_dir = rot_dir / norm( rot_dir ); % renormalize
                % construct matrix
                R_matrix = makehgtform('axisrotate',rot_dir,rot_ang);
                % set hgt matrix
                set(hRotate,'Matrix',MR0*R_matrix);
            end
            %end
        end
        
        % refresh drawing
        drawnow;
    end

% only 1 event on click
    function buttondown_callback(src, ~)
        
        % change status
        mousestatus = 'buttondown';
        % is ctrl pressed?
        ctrl = ismember('control',get(gcbf,'currentmodifier'));
        
        % retrieve the current point
        P = get(gcf,'CurrentPoint');
        % retrieve window geometry
        HGEOM  = get( src, 'Position');
        width  = HGEOM(3);
        height = HGEOM(4);
        
        % SET START POSITION
        if ctrl
            % prepare for translation
            % camera distance
            dCam = norm(get(hAxis,'CameraPosition'));
            % camera view angle
            cva = get(hAxis,'CameraViewAngle');
            % scale from window plane to object plane
            scale = dCam*sind(cva/2)/min([width height]/2);
            START = [scale.*P 0];
            MT0 = get(hTranslate,'Matrix');
            MR0 = get(hRotate,'Matrix');
        else
            % prepare for rotation
            START = point_on_sphere( P, width, height);
            % SET START MATRIX
            MR0 = get(hRotate, 'Matrix');
        end
        % keep ctrl status
        CTRL0 = ctrl;
    end

    function buttonup_callback(~, ~ )
        % change status
        mousestatus = 'buttonup';
        % reset the start position
        START = [0,0,0];
    end

%%%%%%%%%%%% UTILITY FUNCTION %%%%%%%%%%%%%
    function P = point_on_sphere( P, width, height )
        P(3) = 0;
        
        % determine radius of the sphere
        R = min(width, height)/2;
        
        % TRANSFORM the point in window coordinate into
        % the coordinate of a sphere centered in middle of the window
        ORIGIN = [width/2, height/2, 0];
        P = P - ORIGIN;
        
        % normalize position to [-1:1] WRT unit sphere
        % centered at the origin of the window
        P = P / R;
        
        % if position is out of sphere, normalize it to
        % unit length
        L = sqrt(P*P');
        if L > 1
            P = P / L;
            P(3) = 0;
        else
            % add the Z coordinate to the scheme
            P(3) = sqrt( 1 - P(1)^2 - P(2)^2 );
        end
    end

end