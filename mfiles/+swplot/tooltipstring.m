function string = tooltipstring(sObject,swobj)
% generate tooltip string from the data of a graphical object
%
% SWPLOT.TOOLTIPSTRING(sObject,{swobj})
%
% Input:
%
% sObject   Struct, that contains the object data that is clicked on.
% swobj     SpinW object that provides data to the tooltip text, can be
%           empty.
%

if isempty(swobj)
    %       distance energy magnetic-field temperature
    unit = {'' '' '' ''};
else
    unit = swobj.unit.label;
end

% units of energy
unitE = unit{2};
unitR = unit{1};

if ~isempty(unitR)
    unitiR = [unitR '^{-1}'];
else
    unitiR = '';
end


string = '';

switch sObject.name
    case 'atom'
        labelTemp = strword(sObject.label,[1 2],true);
        label1 = labelTemp{1};
        label2 = labelTemp{2};
        string = [label2 ' atom (' label1 ')' newline 'Unit cell:' newline];
        % add cell index and position
        posi = sObject.data(:)';
        string = [string sprintf('[%d, %d, %d]',floor(posi)) newline 'Atomic position:' newline sprintf('[%5.3f, %5.3f, %5.3f]',posi-floor(posi))];
    case 'mag'
        M = sObject.data(:)';
        string = [sprintf('Spin expectation value:\nM = [%5.3f, %5.3f, %5.3f] ',M+1e-6) symbol('hbar')];
    case 'bond'
        M = sObject.data;
        string = [sObject.label ' bond' newline];
        % type      1   Heisenberg exchange
        %           2   Anisotropic exchange
        %           3   DM interaction
        %           4   General matrix
        switch sw_mattype(sObject.data)
            case 1
                % Heisenberg
                string = [string 'Heisenberg exchange' newline 'Value:' newline sprintf('%7.3f ',M(1,1)) ' ' unitE];
            case 2
                % anisotropic exchange
                string = [string 'Anisotropic exchange' newline 'Value:' newline mat2str(diag(M)', [7 3]) ' ' unitE];
            case 3
                % DM interactions
                string = [string 'Antisymmetric exchange (DM vector)' newline 'Value:' newline mat2str([M(2,3) M(3,1) M(1,2)], [7 3]) ' ' unitE];
            case 4
                % General matrix
                string = [string 'General exchange matrix' newline 'Value:' newline mat2str(M, [7 3]) ' ' unitE];
        end
    case {'base' 'base_label'}
        BV = sObject.data;
        RL = 2*pi*inv(BV); %#ok<MINV>
        abc = sqrt(sum(BV.^2,1));
        ang = [sw_angle(BV(:,2),BV(:,3)) sw_angle(BV(:,1),BV(:,3)) sw_angle(BV(:,1),BV(:,2))]*180/pi;
        
        string = ['Lattice parameters' newline...
            sprintf(['%7.3f ' unitR ', '],abc)];
        string = [string(1:(end-2)) newline newline...
            'Angles' newline...
            symbol('alpha') ' = ' sprintf('%7.3f',ang(1)) ' ' symbol('deg') ', '...
            symbol('beta')  ' = ' sprintf('%7.3f',ang(2)) ' ' symbol('deg') ', '...
            symbol('gamma') ' = ' sprintf('%7.3f',ang(3)) ' ' symbol('deg') newline newline...
            'Lattice vectors' newline...
            'a = ' mat2str(BV(:,1)', [7 3]) ' ' unitR newline...
            'b = ' mat2str(BV(:,2)', [7 3]) ' ' unitR newline...
            'c = ' mat2str(BV(:,3)', [7 3]) ' ' unitR newline newline...
            'Reciprocal lattice vectors' newline...
            'h = ' mat2str(RL(1,:), [7 3]) ' ' unitiR newline...
            'k = ' mat2str(RL(2,:), [7 3]) ' ' unitiR newline...
            'l = ' mat2str(RL(3,:), [7 3]) ' ' unitiR];
        
end

end

function string = mat2str(M, format)
% converts matrix into string
%
% MAT2STR(mat, {format}) 
%
% The function converts an NxM matrix into a string with nice zeros, also
% works for matrices with symbolic values.
%
% Input:
%
% M         Two dimensional matrix. format Numerical precision for
%           sprintf function, 1x2 vector: [p1 p2]. It generates the format
%           %p1.p2f. Default is [6 3], optional.
%

if nargin == 1
    p1  = 6;
    p2  = 3;
else
    p1 = format(1);
    p2 = format(2);
end

format = ['%' num2str(p1) '.' num2str(p2) 'f '];
string = [];

for ii = 1:size(M,1)
    string = [string '|']; %#ok<*AGROW>
    for jj = 1:size(M,2)
        if isa(M,'sym')
            string = [string char(M(ii,jj))];
        else
            string = [string sprintf(format,M(ii,jj))];
        end
    end
    if ii < size(M,1)
        string = [string '|' newline];
    else
        string = [string '|'];
    end
end

end