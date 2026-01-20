function varargout = symbolic(obj, symb)
% switches between symbolic/numeric mode
% 
% ### Syntax
% 
% `symb = symbolic(obj)`
%
% `symbolic(obj, symb)`
% 
% ### Description
% 
% `symb = symbolic(obj)` returns `true` if symbolic calculation mode is on,
% `false` for numeric mode.
%  
% `symbolic(obj, symb)` sets whether the calculations are in
% symbolic/numeric (`true`/`false`) mode. Switching to symbolic mode, the
% spin values, matrix elements, magnetic field, magnetic structure and
% physical units are converted into symbolic variables. If this is not
% desired, start with a symbolic mode from the beggining and have full
% control over the values of the above mentioned variables.
% 
% ### See Also
% 
% [spinw] \| [spinw.spinwavesym]
%

% Only returns symb value.
if nargin == 1
    varargout{1} = logical(obj.symb);
    return
end

% No change!
if obj.symb == symb
    return
end

switch symb
    case true
        if ~license('checkout','Symbolic_Toolbox')
            error('spinw:symbolic:NoToolBox','You need Symbolic Math Toolbox installed to run symbolic calculations!');
        end
        
        % Spin values
        obj.unit_cell.S = sym(obj.unit_cell.S);

        % Magnetic structure
        obj.mag_str.F = sym(obj.mag_str.F);
        obj.mag_str.k = sym(obj.mag_str.k);
        
        % Interaction matrices
        nMat = numel(obj.matrix.label);
        if ~isa(obj.matrix.mat,'sym')
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
        
        % Spin values
        symVar1 = symvar(obj.unit_cell.S);
        if ~isempty(symVar1)
            obj.unit_cell.S = double(subs(obj.unit_cell.S,symVar1,ones(1,numel(symVar1))));
        else
            obj.unit_cell.S = double(obj.unit_cell.S);
        end
        
        % Magnetic structure
        symVar1 = symvar(obj.mag_str.F);
        if ~isempty(symVar1)
            obj.mag_str.F = double(subs(obj.mag_str.F,symVar1,ones(1,numel(symVar1))));
        else
            obj.mag_str.F = double(obj.mag_str.F);
        end
        
        symVar1 = symvar(obj.mag_str.k);
        if ~isempty(symVar1)
            obj.mag_str.k = double(subs(obj.mag_str.k,symVar1,ones(1,numel(symVar1))));
        else
            obj.mag_str.k = double(obj.mag_str.k);
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
        
        % Use SI units
        % 0.086173324     Boltzmann constant: k_B [meV/K]
        obj.unit.kB  = 0.086173324;
        % 0.057883818066  Bohr magneton: mu_B [meV/T]
        obj.unit.muB = 0.057883818066;
        
    otherwise
        error('spinw:symbolic:WrongInput','The type of symb input variable has to be logical.')
end

obj.symb = logical(symb);

end