function mouse(hFigure)
% adds mouse callbacks to swplot figure
%
% SWPLOT.MOUSE({hFigure})
%
% Adds rotation and zoom functionality to any swplot figure.
%
% Input:
%
% hFigure       Handle of the swplot figure. Default is the selected
%               figure.
%

if nargin == 0
    hFigure = swplot.activefigure;
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
    
    mousestatus = 'buttonup';
    START = [0 0 0];
    M_previous = get(hRotate,'Matrix');
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
            set(hAxis,'CameraViewAngle',cva*scale);
            set(hAxis,'CameraPosition',get(hAxis,'CameraPosition')/scale);            
        else
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
        ctrl = ismember('control',get(gcbf,'currentmodifier'));
        
        if ctrl
            % pan figure
            hRotate*hTranslate
            
            return
        end
        
        % retrieve the current point
        P = get(gcf,'CurrentPoint');
        % retrieve window geometry
        HGEOM = get(src, 'Position');
        % transform in sphere coordinates (3,4) = (WIDTH, HEIGHT)
        P = point_on_sphere( P, HGEOM(3), HGEOM(4) );
        
        % workaround condition (see point_on_sphere)
        if isnan(P)
            return
        end
        
        %%%%% ARCBALL COMPUTATION %%%%%
        if strcmp(mousestatus, 'buttondown')
            % compute angle and rotation axis
            rot_dir = cross( START, P ); rot_dir = rot_dir / norm( rot_dir );
            rot_ang = acos( dot( P, START ) );
            
            % convert direction in model coordinate system
            rot_dir = M_previous\[rot_dir,0]';
            rot_dir = rot_dir(1:3);
            if norm(rot_dir) > 0
                rot_dir = rot_dir / norm( rot_dir ); % renormalize
                % construct matrix
                R_matrix = makehgtform('axisrotate',rot_dir,rot_ang);
                % set hgt matrix
                set(hRotate,'Matrix',M_previous*R_matrix);
                % refresh drawing
                drawnow;
            end
        end
    end

% only 1 event on click
    function buttondown_callback(src, ~)
        
        % change status
        mousestatus = 'buttondown';
        % retrieve the current point
        P = get(gcf,'CurrentPoint');
        % retrieve window geometry
        HGEOM = get( src, 'Position');
        % SET START POSITION
        START = point_on_sphere( P, HGEOM(3), HGEOM(4) );
        % SET START MATRIX
        M_previous = get(hRotate, 'Matrix');
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
        % the coordinate of a sphere centered in middle window
        ORIGIN = [width/2, height/2, 0];
        P = P - ORIGIN;
        
        % normalize position to [-1:1] WRT unit sphere
        % centered at the origin of the window
        P = P / R;
        
        % if position is out of sphere, normalize it to
        % unit length
        L = sqrt( P*P' );
        if L > 1
            % P = nan; % workaround to stop evaluation
            % disp('out of sphere');
            
            P = P / L;
            P(3) = 0;
        else
            % add the Z coordinate to the scheme
            P(3) = sqrt( 1 - P(1)^2 - P(2)^2 );
        end
    end

end