classdef spinw < handle
    % SPINW class defines data structure and methods to calculate spin wave
    % dispersion in magnetic crystals.
    %
    % obj = SPINW()
    %
    % constructs a new spinw class object, with default parameters.
    %
    % obj = SPINW(obj)
    %
    % constructs new spinw class object. If obj is spinw class, it only
    % checks its data integrity. If obj is struct type, it creates new
    % spinw object and checks data integrity.
    %
    % obj = SPINW(cif_path)
    %
    % construct new spinw class object, where cif_path contains a string of
    % a .cif file path defining an input crystals structure.
    %
    % The data structure behind the spinw object can be accessed by using
    % STRUCT(obj). All fields of the struct type data behind the spinw
    % object are accessible through the main field names of the obj object.
    % For example the lattice parameters:
    %   abc = obj.unit_cell.lat_const;
    %
    % spinw is a handle class, that means that only the handle of the
    % object is copied in a swobj1 = swobj2 command. To create a copy
    % (clone) of an spinw object use:
    %    swobj1 = swobj2.COPY;
    %
    % See also: SPINW.COPY, SPINW.STRUCT,
    % <a href='/Applications/MATLAB_R2012b.app/help/matlab/matlab_oop/comparing-handle-and-value-classes.html'>Comparing handle and value classes</a>.
    %
    % Tutorials and documentation can be found here:
    % <a href='https://www.psi.ch/spinw'>https://www.psi.ch/spinw</a>
    % Forum for questions:
    % <a href='https://groups.google.com/forum/#!forum/spinwforum'>https://groups.google.com/forum/#!forum/spinwforum</a>
    % Lates version and bug reports/feature requests:
    % <a href='https://github.com/tsdev/spinw'>https://github.com/tsdev/spinw</a>
    %
    
    properties (SetObservable)
        % Stores the unit cell parameters.
        % Sub fields are:
        %   'lat_const' lattice constants in a 1x3 vector in Angstrom units
        %   'angle'     (alpha,beta,gamma) angles in 1x3 vector in radian
        %   'sym'       crystal space group, line number in symmetry.dat file
        %
        % See also SPINW.GENLATTICE, SPINW.ABC, SPINW.BASISVECTOR, SPINW.NOSYM.
        lattice
        % Stores the atoms in the crystallographic unit cell.
        % Sub fields are:
        %   'r'         pasitions of the atoms in the unit cell, in a
        %               3 x nAtom matrix, in lattice units
        %   'S'         spin quantum number of the atoms, in a 1 x nAtom
        %               vector, non-magnetic atoms have S=0
        %   'label'     label of the atom, strings in a 1 x nAtom cell
        %   'color'     color of the atom in 3 x nAtom matrix, where every
        %               column is an 0-255 RGB color
        %
        % See also SPINW.ADDATOM, SPINW.ATOM, SPINW.MATOM, SPINW.NEWCELL, SPINW.PLOT.
        unit_cell
        % Stores the crystallographic twin parameters.
        % Sub fields are:
        %   'rotc'      rotation matrices in the xyz coordinate system for
        %               every twin, stored in a 3 x 3 x nTwin matrix
        %   'vol'       volume ratio of the different twins, stored in a
        %               1 x nTwin vector
        %
        % See also SPINW.ADDTWIN, SPINW.TWINQ, SPINW.UNIT_CELL.
        twin
        % Stores 3x3 matrices for using them in the Hailtonian.
        % Sub fields are:
        %   'mat'       stores the actual values of 3x3 matrices, in a
        %               3 x 3 x nMatrix matrix, defult unit is meV
        %   'color'     color assigned for every matrix, stored in a
        %               3 x nMatrix matrix, with 0-255 RGB columns
        %   'label'     label for every matrix, stored as string in a
        %               1 x nMatrix cell
        %
        % See also SPINW.ADDMATRIX, SPINW.NTWIN.
        matrix
        % Stores single ion terms of the Hamiltonian.
        % Sub fields are:
        %   'aniso'     vector contains 1 x nMagAtom integers, each integer
        %               assignes one of the nMatrix from the .matrix field
        %               to a magnetic atom in the spinw.matom list as a single
        %               ion anisotropy (zeros for no anisotropy)
        %   'g'         vector contains 1 x nMagAtom integers, each integer
        %               assignes one of the nMatrix from the .matrix field
        %               to a magnetic atom in the spinw.matom list as a
        %               g-tensor
        %   'field'     external magnetic field stored in a 1x3 vector,
        %               default unit is Tesla
        %   'T'         temperature, scalar, default unit is Kelvin
        %
        % See also SPINW.ADDANISO, SPINW.ADDG, SPINW.GETMATRIX, SPINW.SETMATRIX, SPINW.INTMATRIX.
        single_ion
        % Stores the list of spin-spin couplings.
        % Sub fields are:
        %   'dl'        distance between the unit cells of two interacting
        %               spins, stored in a 3 x nCoupling matrix
        %   'atom1'     first magnetic atom, pointing to the list of
        %               magnetic atoms in spinw.matom list, stored in a
        %               1 x nCoupling vector
        %   'atom2'     second magnetic atom, stored in a  1 x nCoupling
        %               vector
        %   'mat_idx'   stores pointers to matrices for every coupling in a
        %               3 x nCoupling matrix, maximum three matrices per
        %               coupling (zeros for no coupling)
        %   'idx'       increasing indices for the symmetry equivalent
        %               couplings, starting with 1,2,3...
        %   'type'      Type of coupling corresponding to mat_idx matrices.
        %               Default is 0 for quadratic exchange. type = 1 for
        %               biquadratic exchange.
        %
        % See also SPINW.GENCOUPLING, SPINW.ADDCOUPLING, SPINW.FIELD.
        coupling
        % Stores the magnetic structure.
        % Sub fields are:
        %   'S'         stores the moment direction for every spin in the
        %               crystallographic or magnetic supercell in a
        %               3 x nMagExt matrix, where nMagExt = nMagAtom*prod(N_ext)
        %   'k'         magnetic ordering wave vector in a 3x1 vector
        %   'n'         normal vector to the rotation of the moments in
        %               case of non-zero ordering wave vector, dimensions
        %               are 3x1
        %   'N_ext'     Size of the magnetic supercell, default is [1 1 1]
        %               if the magnetic cell is identical to the
        %               crystallographic cell, the 1x3 vector extends the
        %               cell along the a, b and c axis
        %   'rdip'      Maximum distance until the dipolar interaction is
        %               calculated. Zero vlue means no dipolar interactions
        %               are considered.
        %
        % See also SPINW.GENMAGSTR, SPINW.OPTMAGSTR, SPINW.ANNEAL, SPINW.MOMENT, SPINW.NMAGEXT, SPINW.STRUCTFACT.
        mag_str
        % Stores the physical units in the Hamiltonian.
        % Defaults are meV, Tesla Angstrom and Kelvin.
        % Sub fields are:
        %   'kB'        Boltzmann constant, default is 0.0862 [meV/K]
        %   'muB'       Bohr magneton, default is 0.0579 [meV/T]
        %   'mu0'       Vacuum permeability, 201.335431 [T^2*Angstrom^3/meV]
        unit
        % Stores the cache, it should be only used to check consistency of
        % the code. The stored values should not be changed by the user!
        % Sub fields are:
        %   'matom'     Data on the magnetic unit cell.
        cache = struct('matom',[]);
        
    end
    
    properties (Access = private)
        propl         % stores the property change listener handles
        sym  = false; % stores whether the couplings are generated under symmetry constraints
        symb = false; % stores whether the calculation are done symbolically
        fid  = 1;     % stores the file ID of the text output, default is the Command Window
        ver    = sw_version;
    end
    
    methods
        function obj = spinw(varargin)
            % SPINW constructor
            %
            
            if nargin==0
                objS = initfield(struct);
                fNames = fieldnames(objS);
                for ii = 1:length(fNames)
                    obj.(fNames{ii}) = objS.(fNames{ii});
                end
                return;
            end
            
            firstArg = varargin{1};
            if isa(firstArg, 'spinw') %  It is used when objects are passed as arguments.
                obj = copy(firstArg);
                return;
            end
            
            if isstruct(firstArg)
                objS = initfield(firstArg);
                
                % change lattice object
                if size(objS.lattice.lat_const,1)==3
                    objS.lattice.lat_const = objS.lattice.lat_const';
                end
                if size(objS.lattice.angle,1)==3
                    objS.lattice.angle = objS.lattice.angle';
                end
                
                validate(objS);
                fNames = fieldnames(objS);
                for ii = 1:length(fNames)
                    obj.(fNames{ii}) = objS.(fNames{ii});
                end
                return;
            end
            if ischar(firstArg)
                % import data from cif file
                
                objS = initfield(struct);
                
                cif0 = cif(firstArg);
                fprintf0(1,'Crystal structure is imported from %s.\n',firstArg);
                abc0 = [cif0.cell_length_a cif0.cell_length_b cif0.cell_length_c];
                ang0 = [cif0.cell_angle_alpha cif0.cell_angle_beta cif0.cell_angle_gamma];
                sym0 = cif0.('symmetry_space_group_name_H-M');
                %symi0 = cif0.symmetry_Int_Tables_number;
                
                if ~isempty(cif0.symmetry_equiv_pos_as_xyz)
                    xyz0 = cif0.symmetry_equiv_pos_as_xyz;
                elseif ~isempty(cif0.space_group_symop_operation_xyz)
                    xyz0 = cif0.space_group_symop_operation_xyz;
                else
                    warning('spinw:WrongFormat','Missing symmetry operators, using P1')
                    xyz0 = {'x,y,z'};
                end
                
                xyz0 = sprintf('%s; ',xyz0{:}); xyz0 = xyz0(1:end-2);
                %name0 = cif0.atom_site_type_symbol';
                cell0=[cif0.atom_site_label cif0.atom_site_type_symbol];
                name0 = cellfun(@(x,y)strjoin({x y}),cell0(:,1),cell0(:,2),'UniformOutput',false)';
                r0 = mod([cif0.atom_site_fract_x cif0.atom_site_fract_y cif0.atom_site_fract_z]',1);
                %occ0 = cif0.atom_site_occupancy;
                
                if numel(abc0)==3
                    objS.lattice.lat_const = abc0;
                end
                if numel(ang0) == 3
                    objS.lattice.angle = ang0*pi/180;
                end
                if numel(xyz0) > 3
                    symIdx = sw_addsym(xyz0,sym0);
                    objS.lattice.sym = int32(symIdx);
                    
                end
                if size(name0,2) == size(r0,2)
                    nAtom = numel(name0);
                    
                    objS.unit_cell.r = r0;
                    objS.unit_cell.label = name0;
                    
                    col = zeros(3,nAtom);
                    for ii = 1:nAtom
                        aName = strword(name0{ii},2,true);
                        col0 = sw_atomdata(aName{1}(aName{1}>57),'color')';
                        col(:,ii) = col0;
                        [~, ~, objS.unit_cell.S(ii)] = sw_mff(aName{1});
                    end
                    objS.unit_cell.color = int32(col);
                end
                
                fNames = fieldnames(objS);
                for ii = 1:length(fNames)
                    obj.(fNames{ii}) = objS.(fNames{ii});
                end
                
            end
            
            
        end % .spinw
        
        function nMagExt = nmagext(obj)
            % gives the number of magnetic atoms in the magnetic supercell
            %
            % nMagExt = NMAGEXT(obj)
            %
            
            nMagExt = size(obj.mag_str.S,2);
        end
        function nAtom = natom(obj)
            % gives the number of symmetry unrelated atoms in the unit cell
            %
            % nAtom = NATOM(obj)
            %
            
            nAtom = size(obj.unit_cell.r,2);
        end
        function nBond = nbond(obj)
            % gives the number of bonds defined in the spinw object
            %
            % nBond = NBOND(obj)
            %
            
            nBond = size(obj.coupling.idx,2);
        end
        
        function nMat = nmat(obj)
            % gives the number of matrices defined in an spinw object
            %
            % nMat = NMAT(obj)
            %
            
            nMat = size(obj.matrix.mat,3);
        end
        function nTwin = ntwin(obj)
            % gives the number of twins
            %
            % nTwin = NTWIN(obj)
            %
            
            nTwin = size(obj.twin.vol,2);
        end
                
    end
    
    methods(Hidden=true)
        function modmatom(obj, ~, ~)
            % listening to the change of the lattice or unit_cell fields
            
            % delete the stored magnetic atom positions
            obj.cache.matom = [];
            % remove the listener
            delete(obj.propl);
            % fprintf('Property changed!\n')
        end
        function lh = addlistener(varargin)
            lh = addlistener@handle(varargin{:});
        end
        function notify(varargin)
            notify@handle(varargin{:});
        end
        function Hmatch = findobj(varargin)
            Hmatch = findobj@handle(varargin{:});
        end
        function p = findprop(varargin)
            p = findprop@handle(varargin{:});
        end
        function TF = eq(varargin)
            TF = eq@handle(varargin{:});
        end
        function TF = ne(varargin)
            TF = ne@handle(varargin{:});
        end
        function TF = lt(varargin)
            TF = lt@handle(varargin{:});
        end
        function TF = le(varargin)
            TF = le@handle(varargin{:});
        end
        function TF = gt(varargin)
            TF = gt@handle(varargin{:});
        end
        function TF = ge(varargin)
            TF = ge@handle(varargin{:});
        end
    end % classdef
    
end