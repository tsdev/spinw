function view(ax,hFigure)
% control 3D view of swplot
%
% SWPLOT.VIEW(ax,{hFigure})
%
% Input:
%
% ax        String to change view of the 3D plot. Possible options:
%               'ab','bc','ac'  the two axes define the view plane,
%               'hk','kl','hl'  the two reciprocal lattice vectors define
%                               the view plane.
% hFigure   Handle of the swplot figure window, optional.
%

if nargin<2
    % find active figure
    hFigure = swplot.activefigure;
end

% hgtransform object
h = getappdata(hFigure,'h');

if isempty(h)
    return
end

% matrix that does the rotation
%hRot = get(h,'Parent');
% set mouse rotation to zero
%set(get(h,'Parent'),'Matrix',eye(4));

bv = getappdata(hFigure,'base');

if any(ax>'c')
    % for hkl
    bv = inv(bv)'*2*pi;
end

bvN = bsxfun(@rdivide,bv,sqrt(sum(bv.^2,1)));

x = bvN(:,1);
y = bvN(:,2);
z = bvN(:,3);

cc(1) = y'*z;
cc(2) = x'*z;
cc(3) = x'*y;

switch ax
    case {'ab' 'c' 'hk' 'l'}
        % ab
        xp = [1;0;0];
        yp = [cc(3);sqrt(1-cc(3)^2);0];
        zp = [cc(2);(cc(1)-yp(1)*cc(2))/yp(2)];
        zp(3) = sqrt(1-sum(zp.^2));
        zp(3) = sign(cross(xp,yp)'*zp)*zp(3);
    case {'bc' 'a' 'kl' 'h'}
        % bc
        yp = [1;0;0];
        zp = [cc(1);sqrt(1-cc(1)^2);0];
        xp = [cc(3);(cc(2)-zp(1)*cc(3))/zp(2)];
        xp(3) = sqrt(1-sum(xp.^2));
        xp(3) = sign(cross(xp,yp)'*zp)*xp(3);
    case {'ac' 'b' 'hl' 'k'}
        % ac
        xp = [1;0;0];
        zp = [cc(2);sqrt(1-cc(2)^2);0];
        yp = [cc(3);(cc(1)-zp(1)*cc(3))/zp(2)];
        yp(3) = sqrt(1-sum(yp.^2));
        yp(3) = sign(cross(xp,yp)'*zp)*yp(3);
    otherwise
        error('view:WrongInput','Wrong input!')
end

M = swplot.transform(hFigure);
M(1:3,1:3) = [xp yp zp]/bvN;
swplot.transform(M,hFigure);

end