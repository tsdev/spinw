function varargout = symbolic(obj, symb)
% change between symbolic/numeric calculation
%
% symb = SYMBOLIC(obj)
%
% Returns true is symbolic calculation mode is on, false for numeric mode.
%
% SYMBOLIC(obj, symb)
%
% symb sets whether the calculations are symbolic/numeric (true/false).
%
% See also SW, SW.SPINWAVESYM.
%

% Only returns symb value.
if nargin == 1
    varargout{1} = logical(obj.symb);
    return;
end

% No change!
if obj.symb == symb
    return
end

switch symb
    case true
        v = ver;
        if ~any(strcmp('Symbolic Math Toolbox', {v.Name}))
            error('sw:symbolic:NoToolBox','You need Symbolic Math Toolbox installed to run symbolic calculations!');
        end
        
        % Magnetic structure
        obj.mag_str.S = sym(obj.mag_str.S);
        obj.mag_str.k = sym(obj.mag_str.k);
        obj.mag_str.n = sym(obj.mag_str.n);
        
        % Interaction matrices
        nMat = numel(obj.matrix.label);
        if isa(obj.matrix.mat,'sym')
            % matrices are already symbolic
            %elseif nMat == 0
            %    obj.matrix.mat = sym(obj.matrix.mat);
            %else
            mat0 = sym(obj.matrix.mat*0);
            for ii = 1:nMat
                symVar = sym(obj.matrix.label{ii},'real');
                mat0(:,:,ii) = obj.matrix.mat(:,:,ii)*symVar;
            end
            obj.matrix.mat = mat0;
        end
        
        % Magnetic field
        if ~isa(obj.field,'sym')
            obj.single_ion.field = obj.single_ion.field * sym('B','real');
        end
        
        % Units
        obj.unit.kB  = sym('kB','positive');
        obj.unit.muB = sym('muB','positive');
        
    case false
        % Create double type properties
        
        % Magnetic structure
        symVar1 = symvar(obj.mag_str.S);
        if ~isempty(symVar1)
            obj.mag_str.S = double(subs(obj.mag_str.S,symVar1,ones(1,numel(symVar1))));
        else
            obj.mag_str.S = double(obj.mag_str.S);
        end
        
        symVar1 = symvar(obj.mag_str.k);
        if ~isempty(symVar1)
            obj.mag_str.k = double(subs(obj.mag_str.k,symVar1,ones(1,numel(symVar1))));
        else
            obj.mag_str.k = double(obj.mag_str.k);
        end
        
        symVar1 = symvar(obj.mag_str.n);
        if ~isempty(symVar1)
            obj.mag_str.n = double(subs(obj.mag_str.n,symVar1,ones(1,numel(symVar1))));
        else
            obj.mag_str.n = double(obj.mag_str.n);
        end
        
        % Interaction matrices
        symVar1 = symvar(obj.matrix.mat);
        if ~isempty(symVar1)
            obj.matrix.mat = double(subs(obj.matrix.mat,symVar1,ones(1,numel(symVar1))));
        else
            obj.matrix.mat = double(obj.matrix.mat);
        end
        
        % Magnetic field
        symVar1 = symvar(obj.single_ion.field);
        if ~isempty(symVar1)
            obj.single_ion.field = double(subs(obj.single_ion.field,symVar1,ones(1,numel(symVar1))));
        else
            obj.single_ion.field = double(obj.single_ion.field);
        end
        
        % Units
        obj.unit.kB  = 0.086173324;
        obj.unit.muB = 0.057883818066;
        
    otherwise
        error('sw:symbolic:WrongInput','The type of symb input variable has to be logical.')
end

obj.symb = logical(symb);

end