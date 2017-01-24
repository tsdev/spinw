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
        a2Idx = sObject.data(end);
        aIdx  = swobj.atom.idx(a2Idx);
        maIdx = sum(swobj.atom.mag(1:a2Idx));
        
        labelTemp = strword(sObject.label,[1 2],true);
        %label1 = labelTemp{1};
        label2 = labelTemp{2};
        string = [label2 ' atom (''' sObject.label ''')' newline];
        if swobj.atom.mag(a2Idx)
            string = [string sprintf('Magnetic, S = %g',swobj.unit_cell.S(aIdx)) newline];
        else
           string = [string 'Non-magnetic' newline];
        end
        string = [string sprintf('Index in the spinw.unit\\_cell list: #%d',aIdx) newline];
        string = [string sprintf('Index in the spinw.atom() list: #%d',a2Idx) newline];
        if swobj.atom.mag(a2Idx)
            string = [string sprintf('Index in the spinw.matom() list: #%d',maIdx) newline];
        end
        
        string = [string 'Unit cell:' newline];
        % add cell index and position
        cellindex = sObject.data(1:3)';
        pos = swobj.unit_cell.r(:,aIdx);
        
        string = [string sprintf('[%d, %d, %d]',cellindex) newline 'Atomic position:' newline sprintf('[%5.3f, %5.3f, %5.3f]',pos)];
    case 'mag'
        Mplot = sObject.data(1:3)';
        pos   = sObject.data(4:6)';
        cellindex = floor(pos);
        maIdx = sObject.data(7);
        
        % position in the magnetic list
        nExt = double(swobj.mag_str.nExt);
        nMagAtom = numel(swobj.matom.idx);
        
        mIdx = nMagAtom*(cellindex(1) + cellindex(2)*nExt(1) + cellindex(3)*prod(nExt(1:2)))+maIdx;
        
        string = ['Magnetic moment vector' newline];
       
        string = [string sprintf('Size of magnetic supercell: [%d,%d,%d]',swobj.mag_str.nExt) newline];
        if swobj.magstr.exact
            string = [string 'The supercell is ideal: S=F!' newline];
        else
            string = [string 'The supercell forces an approximation: S~F' newline];
        end
        
        nK = size(swobj.mag_str.k,2);
        if nK == 1
            string = [string sprintf('Magnetic propagation vector:\n[%5.3f, %5.3f, %5.3f]',swobj.mag_str.k) newline];
        else
            string = [string 'Magnetic propagation vectors:' newline];
            for ii = 1:nK
                string = [string sprintf('[%5.3f, %5.3f, %5.3f]\n',swobj.mag_str.k(:,ii))];
            end
        end
        string = [string sprintf('Moment index in spinw.mag\\_str: #%d',mIdx) newline];
        string = [string sprintf('Spin expectation value:\nM = [%5.3f, %5.3f, %5.3f] ',Mplot+1e-6) symbol('hbar') newline];
        string = [string sprintf('Normalized value:\nM = [%5.3f, %5.3f, %5.3f] ',(Mplot+1e-6)/norm(Mplot))];
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