function ellM = sw_quadell(M, toplot)
% calculates and plots the parameters of an ellipsoid from a quadratic form
%
% ellM = SW_QUADELL(M, plot)
%
% It calculates the parameters of an ellipsoid belonging to the input
% quadratic form and plots it. Can be used as a tool to visualise
% anisotropy matrices, g-tensors, etc.
%
% Input:
%
% M         Quadratic form coefficients in a 3x3 matrix.
% toplot    If true, the ellipsoid is plotted. Default is true.
%
% Output:
%
% ellM      Contais the principal axis of the ellipsoid in the columns, 3x3
%           matrix.
%
% On the plot the principal axes are shown by a red line, the x,y and z
% axes are shown by black line.
%

if nargin == 0
    help sw_quadell;
    return
end

if nargin < 2
    toplot = true;
end

% take only the symmetric part of M
if any(any(M-M'))
    M = (M+M')/2;
    warning('sw_quadell:Antisym','M is not symmetric, throwing away the antisymmetric components.');
end

% minimum ellipsoid size
epsilon = 0.1;
% maximum radius of the ellipsoid
ellAniso = 1;
% center
c = [0 0 0];
% number of surface points
N = 50;


% Calculates the main radiuses of the ellipsoid.
[V, R] = eig(M);

% Creates positive definite matrix by adding constant to all
% eigenvalues.
R0      = diag(R);
epsilon = sqrt(sum(R0.^2))*epsilon;
%dR      = 1./(R0-min(R0)+epsilon);
dR      = (R0-min(R0)+epsilon);
dR0     = max(dR);
if dR0 ~= 0
    dR = dR/dR0;
else
    dR = [1 1 1];
end

R = diag(dR)*ellAniso;

ell.mat = V*R;

ellM = ell.mat;

if toplot
    figure;
    [spx, spy, spz] = sphere(N);
    
    ell.xyz = ell.mat*[spx(:) spy(:) spz(:)]';
    % draw main circles
    c1 = ell.mat*sw_circle([0 0 0]',[1 0 0]',1.01,N);
    c2 = ell.mat*sw_circle([0 0 0]',[0 1 0]',1.01,N);
    c3 = ell.mat*sw_circle([0 0 0]',[0 0 1]',1.01,N);
    
    plot3(c1(1,:),c1(2,:),c1(3,:),'color','k','linewidth',3);
    hold on
    plot3(c2(1,:),c2(2,:),c2(3,:),'color','k','linewidth',3);
    plot3(c3(1,:),c3(2,:),c3(3,:),'color','k','linewidth',3);
    
    ell.x = reshape(ell.xyz(1,:),[1 1]*N+1);
    ell.y = reshape(ell.xyz(2,:),[1 1]*N+1);
    ell.z = reshape(ell.xyz(3,:),[1 1]*N+1);
    
    ell.h = surf(ell.x+c(1),ell.y+c(2),ell.z+c(3));
    set(ell.h,'EdgeAlpha',0);
    
    set(gca,'Position',[0 0 1 1]);
    set(gca,'Color','none');
    set(gca,'Box','off');
    set(gca,'Clipping','Off');
    daspect([1 1 1]);
    pbaspect([1 1 1]);
    axis on
    axis vis3d
    material dull;
    grid on
    xlabel('x')
    ylabel('y')
    zlabel('z')
    
    aLim = axis;
    line(aLim(1:2),[0 0],[0 0],'color','k','linewidth',2);
    line([0 0],aLim(3:4),[0 0],'color','k','linewidth',2);
    line([0 0],[0 0],aLim(5:6),'color','k','linewidth',2);
    
    aLine = ell.mat(:,1); aLine = aLine/norm(aLine);
    bLine = ell.mat(:,2); bLine = bLine/norm(bLine);
    cLine = ell.mat(:,3); cLine = cLine/norm(cLine);
    
    lLim = [-1 1]*max(max(abs(aLim)));
    line(aLine(1)*lLim,aLine(2)*lLim,aLine(3)*lLim,'color','r','linewidth',4)
    line(bLine(1)*lLim,bLine(2)*lLim,bLine(3)*lLim,'color','r','linewidth',4)
    line(cLine(1)*lLim,cLine(2)*lLim,cLine(3)*lLim,'color','r','linewidth',4)
    colormap([1 1 1]/2)
end

end

function points = sw_circle(r0, n, R, N)
% creates the 3D coordinates of the circle circumference
% 
% ### Syntax
% 
% `R = sw_circle(r0, n, r, n)`
% 
% ### Description
%
% `points = sw_circle(r0, n, r, n)` generates the 3D coordinates of a
% circle circumference defined by the position of the circle, normal vector
% and radius.
%
% ### Input Arguments
%
% `r0`
% : Center of circle in a column vector with 3 elements.
%
% `n`
% : Normal to the circle surface, in a column vector with 3 elements.
%
% `R`
% : Radius of the circle.
% 
% `N`
% : Number of points on the circumference.
% 
% ### Output Arguments
%
% `R`
% : Matrix with dimensions of $[3\times N]$ containing the 3D coordinates
%   in columns.
%

if nargin == 0
    help sw_circle
    return
end

if any(cross(n,[0; 0; 1]))
    a = cross(n,[0; 0; 1]);
else
    a = cross(n,[0; 1; 0]);
end

b = cross(n,a);
a = a/norm(a);
b = b/norm(b);

phi = linspace(0,2*pi,N);

points = R*(a*cos(phi)+b*sin(phi))+repmat(r0,1,N);

end