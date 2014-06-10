classdef (ConstructOnLoad) sw < class_handlelight
    % SW class defines data structure and methods to calculate spin wave
    % dispersion in magnetic crystals.
    %
    % obj = SW() 
    %
    % constructs a new sw class object, with default parameters.
    %
    % obj = SW(obj)
    %
    % constructs new sw class object. If obj is sw class, it only checks
    % its data integrity. If obj is struct type, it creates new sw object
    % and checks data integrity.
    %
    % obj = SW(cif_path) 
    %
    % construct new sw class object, where cif_path contains a string of a
    % .cif file path defining an input crystals structure.
    %
    % The data structure behind the sw object can be accessed by using
    % STRUCT(obj). All fields of the struct type data behind the sw object
    % are accessible through the main field names of the obj object. For
    % example the lattice parameters:
    %   abc = obj.unit_cell.lat_const;
    %
    % sw is a handle class, that means that only the handle of the object
    % is copied in a swobj1 = swobj2 command. To create a copy (clone) of
    % an sw object use:
    %    swobj1 = swobj2.COPY;
    %
    % See also: SW.COPY, SW.STRUCT,
    % <a href='/Applications/MATLAB_R2012b.app/help/matlab/matlab_oop/comparing-handle-and-value-classes.html'>Comparing handle and value classes</a>.
    %
    % Information in a blog form:
    % <a href='http://spinw.tumblr.com'>http://spinw.tumblr.com</a>
    % Documentation can be found here:
    % <a href='https://wiki.helmholtz-berlin.de/spinw'>https://wiki.helmholtz-berlin.de/spinw</a>
    % Forum for questions:
    % <a href='https://groups.google.com/forum/#!forum/spinwforum'>https://groups.google.com/forum/#!forum/spinwforum</a>
    % Lates version and bug reports/feature requests:
    % <a href='http://code.google.com/p/spinw/'>http://code.google.com/p/spinw/</a>
    %
    
    properties
        % Field stores the crystallographic unit cell parameters.
        % Sub fields are:
        %   'lat_const' lattice constants in a 1x3 vector in Angstrom units
        %   'angle'     (alpha,beta,gamma) angles in 1x3 vector in radian
        %   'sym'       crystal space group, line number in symmetry.dat file
        %
        % See also SW.GENLATTICE, SW.ABC, SW.BASISVECTOR, SW.NOSYM.
        lattice
        % Field stores the atoms in the crystallographic unit cell.
        % Sub fields are:
        %   'r'         pasitions of the atoms in the unit cell, in a
        %               3 x nAtom matrix, in lattice units
        %   'S'         spin quantum number of the atoms, in a 1 x nAtom
        %               vector, non-magnetic atoms have S=0
        %   'label'     label of the atom, strings in a 1 x nAtom cell
        %   'color'     color of the atom in 3 x nAtom matrix, where every
        %               column is an 0-255 RGB color
        %
        % See also SW.ADDATOM, SW.ATOM, SW.MATOM, SW.NEWCELL, SW.PLOT.
        unit_cell
        % Field stores crystallographic twins.
        % Sub fields are:
        %   'rotc'      rotation matrices in the xyz coordinate system for
        %               every twin, stored in a 3 x 3 x nTwin matrix
        %   'vol'       volume ratio of the different twins, stored in a
        %               1 x nTwin vector
        %
        % See also SW.ADDTWIN, SW.TWINQ, SW.UNIT_CELL.
        twin
        % Field stores 3x3 matrices for using them in the Hailtonian.
        % Sub fields are:
        %   'mat'       stores the actual values of 3x3 matrices, in a
        %               3 x 3 x nMatrix matrix, defult unit is meV
        %   'color'     color assigned for every matrix, stored in a
        %               3 x nMatrix matrix, with 0-255 RGB columns
        %   'label'     label for every matrix, stored as string in a
        %               1 x nMatrix cell
        %
        % See also SW.ADDMATRIX, SW.NTWIN.
        matrix
        % Field stores single ion terms of the Hamiltonian.
        % Sub fields are:
        %   'aniso'     vector contains 1 x nMagAtom integers, each integer
        %               assignes one of the nMatrix from the .matrix field
        %               to a magnetic atom in the sw.matom list as a single
        %               ion anisotropy (zeros for no anisotropy)
        %   'g'         vector contains 1 x nMagAtom integers, each integer
        %               assignes one of the nMatrix from the .matrix field
        %               to a magnetic atom in the sw.matom list as a
        %               g-tensor
        %   'field'     external magnetic field stored in a 1x3 vector,
        %               default unit is Tesla
        %   'T'         temperature, scalar, default unit is Kelvin
        %
        % See also SW.ADDANISO, SW.ADDG, SW.GETMATRIX, SW.SETMATRIX, SW.INTMATRIX.
        single_ion
        % Field stores the list of spin-spin couplings.
        % Sub fields are:
        %   'dl'        distance between the unit cells of two interacting
        %               spins, stored in a 3 x nCoupling matrix
        %   'atom1'     first magnetic atom, pointing to the list of
        %               magnetic atoms in sw.matom list, stored in a
        %               1 x nCoupling vector
        %   'atom2'     second magnetic atom, stored in a  1 x nCoupling
        %               vector
        %   'mat_idx'   stores pointers to matrices for every coupling in a
        %               3 x nCoupling matrix, maximum three matrices per
        %               coupling (zeros for no coupling)
        %   'idx'       increasing indices for the symmetry equivalent
        %               couplings, starting with 1,2,3...
        %
        % See also SW.GENCOUPLING, SW.ADDCOUPLING, SW.FIELD.
        coupling
        % Field stores the magnetic structure.
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
        %
        % See also SW.GENMAGSTR, SW.OPTMAGSTR, SW.ANNEAL, SW.MOMENT, SW.NMAGEXT, SW.STRUCTFACT.
        mag_str
        % Field stores the physical units in the Hamiltonian. Default are
        % meV, Tesla and Kelvin.
        % Sub fields are:
        %   'kB'        Boltzmann constant, default is 0.0862
        %   'muB'       Bohr magneton, default is 0.0579
        unit
    end
    
    properties (Access = private)
        matomstore = [];
        sym  = true;  % stores whether the couplings are generated under symmetry constraints
        symb = false; % stores whether the calculation are done symbolically
        fid  = 1;     % stores the file ID of the text output, default is the Command Window
        Elabel = 'meV';
        Qlabel = 'Angstrom^{-1}';
        Rlabel = 'Angstrom';
        Blabel = 'T';
        Tlabel = 'K';
    end
    
    methods
        function obj = sw(varargin)
            % SW constructor
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
            if isa(firstArg, 'sw') %  It is used when objects are passed as arguments.
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
                xyz0 = cif0.symmetry_equiv_pos_as_xyz; xyz0 = sprintf('%s; ',xyz0{:}); xyz0 = xyz0(1:end-2);
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
                        [~, ~, objS.unit_cell.S(ii)] = sw_mff(name0{ii});
                    end
                    objS.unit_cell.color = int32(col);
                end
                
                fNames = fieldnames(objS);
                for ii = 1:length(fNames)
                    obj.(fNames{ii}) = objS.(fNames{ii});
                end
                
            end
            
            
        end % .sw
        
        function objC = copy(obj)
            % clones sw object with all data
            %
            % newObj = COPY(obj)
            %
            % Use this command instead of the '=' sign if you want to
            % create an independent duplicate of an sw class object.
            %
            % Input:
            %
            % obj       sw class object.
            %
            % Output:
            %
            % newObj    New sw class object that contains all the data of
            %           obj.
            %
            % Example:
            %
            % cryst = sw;
            % cryst.addmatrix('label','J1','value',3.1415)
            %
            % cryst1 = cryst;
            % cryst2 = cryst.copy;
            %
            % cryst.addmatrix('label','J1','value',1)
            % J1a = cryst1.matrix.mat;
            % J1b = cryst2.matrix.mat;
            %
            % Where J1a will be a matrix with 1 in the diagonal, while J1b
            % has 3.1415 in the diagonal. If cryst is changed, cryst1 will
            % be changed as well and viece versa, since they point to the
            % same object. However cryst2 is independent of cryst.
            %
            % See also SW, SW.STRUCT.
            %
            
            objS = struct(obj);
            objC = sw(objS);
            
            % copy the private properties
            objC.sym    = obj.sym;
            objC.symb   = obj.symb;
            objC.fid    = obj.fid;
            objC.Elabel = obj.Elabel;
            objC.Qlabel = obj.Qlabel;
            objC.Rlabel = obj.Rlabel;
            objC.Blabel = obj.Blabel;
            objC.Tlabel = obj.Tlabel;
            
        end % copy

        function abc = abc(obj)
            % returns lattice parameters and angles
            %
            % latVect = ABC(obj)
            %
            % Input:
            %
            % obj       sw class object
            %
            % Output:
            %
            % latVect   Vetor with elements [a, b, c, alpha, beta, gamma],
            %           contains the lattice parameters and angles in
            %           Angstrom and degree units respectively.
            %
            % See also SW.HORACE.
            %
            
            abc = [obj.lattice.lat_const obj.lattice.angle*180/pi];
        end
        function nMagExt = nmagext(obj)
            % gives the number of magnetic atoms in the magnetic supercell
            %
            % nMagExt = NMAGEXT(obj)
            %
            
            nMagExt = size(obj.mag_str.S,2);
        end
        function nTwin = ntwin(obj)
            % gives the number of twins
            %
            % nTwin = NTWIN(obj)
            %
            
            nTwin = size(obj.twin.vol,2);
        end
        
        function varargout = temperature(obj,varargin)
            % get/set stored temperature value
            %
            % TEMPERATURE(obj, T)
            %
            % If T is defined, it sets the temperature stored in obj object
            % to T, where T is scalar. The units of temerature is
            % determined by the sw.unit.kB value, default is Kelvin.
            %
            % T = TEMPERATURE(obj)
            %
            % The function returns the current temperature value stored in
            % obj.
            %
            
            if nargin == 1
                varargout{1} = obj.single_ion.T;
            elseif nargin == 2
                T = varargin{1};
                if numel(T) == 1
                    obj.single_ion.T = T;
                else
                    error('sw:temperature:ArraySize','Input temperature has to be scalar!');
                end
                if nargout > 0
                    varargout{1} = obj;
                end
            end
            
        end % .temperature
        
        function sym = symmetry(obj)
            % true if space group is used to generate couplings
            %
            % sym = SYMMETRY(obj)
            %
            % If true, equivalent couplings are generated based on the
            % crystal space group and all matrices (interaction, anisotropy
            % and g-tensor) are transformed according to the symmetry
            % operators. If false, equivalent couplings are generated based
            % on bond length, equivalent matrices won't be transformed
            % (all identical).
            %
            % To change it use sw.gencoupling with the forceNoSym option.
            % To remove all symmetry operators use sw.nosym.
            %
            % See also SW, SW.NOSYM, SW.GENCOUPLING.
            %
            
            sym = obj.sym;
            
        end % .symmetry

        function fidOut = fileid(obj,fid)
            % determines where the text out is written
            %
            % FILEID(obj,fid)
            %
            % Determines the text output of all sw class methods. Default
            % is 1, where all output is printed onto the MATLAB Command
            % Window.
            %
            % fidOut = FILEID(obj)
            %
            % Outputs the stored fileID value.
            %
            % See also SW, SW.DISPLAY.
            %
            if nargin > 1
                obj.fid = fid;
            end
            
            if nargout > 0
                fidOut = obj.fid;
            end
            
        end % .fileid

        function varargout = notwin(obj)
            % removes any twin added to the sw object
            %
            % NOTWIN(obj)
            %
            % The function keeps only the original twin.
            %
            
            obj.twin.vol = 1;
            obj.twin.rotc = eye(3);
            
            if nargout >0
                varargout{1} = obj;
            end
        end
            
        function varargout = symbolic(obj, symb)
            % true/false for symbolic/numerical calculation
            %
            % symb = SYMBOLIC(obj)
            %
            % If true, magnetic structure are spin wave dispersions are
            % calculated symbolically.
            %
            % SYMBOLIC(obj, symb)
            %
            % symb sets whether the calculations are symbolic/numerical
            % (true/false).
            %
            % See also SW, SW.SPINWAVESYM.
            %
            
            if nargin == 1
                symb = obj.symb;
                
                if symb
                    v = ver;
                    if ~any(strcmp('Symbolic Math Toolbox', {v.Name}))
                        symb = false;
                        warning('You need Symbolic Math Toolbox installed to run symbolic calculations!')
                    end
                end
                
            elseif symb == true
                
                obj.mag_str.S = sym(obj.mag_str.S); %#ok<*CPROP>
                obj.mag_str.k = sym(obj.mag_str.k);
                obj.mag_str.n = sym(obj.mag_str.n);
                
                nMat = numel(obj.matrix.label);
                if isa(obj.matrix.mat,'sym')
                    % matrices are already symbolic
                elseif nMat == 0
                    obj.matrix.mat = sym(obj.matrix.mat);
                else
                    mat0 = sym(obj.matrix.mat*0);
                    for ii = 1:nMat
                        symVar = sym(obj.matrix.label{ii},'real');
                        mat0(:,:,ii) = obj.matrix.mat(:,:,ii)*symVar;
                    end
                    obj.matrix.mat = mat0;
                end
                
                obj.symb = true;
            end
            
            if (nargin < 2) || (nargout > 0)
                varargout = {logical(symb)};
            end
            
        end % .symbolic
        
    end
    methods (Hidden=true)
        function obj = addmagfield(obj, B)
            warning('sw:addmagfield:Obsolete','addmagfield is obsolete, use sw.field([Bx By Bz]) instead!');
            obj.single_ion.field = B;
        end % .addmagfield
        function obj = sw_gencoupling(obj, varargin)
            warning('sw:sw_gencoupling:Obsolete','sw_gencoupling is obsolete, use gencoupling instead!')
            obj = gencoupling(obj, varargin{:});
        end % sw_gencoupling
        
        function obj = addjtype(obj, varargin)
            warning('sw:addjtype:Obsolete','addjtype is obsolete, use addmatrix instead!')
            obj = addmatrix(obj, varargin{:});
        end
        function obj = addj(obj, varargin)
            warning('sw:addj:Obsolete','addj is obsolete, use addcoupling instead!')
            obj = addcoupling(obj, varargin{:});
        end
        function obj = sw_magstr(obj, varargin)
            warning('sw:sw_magstr:Obsolete','sw_magstr is obsolete, use genmagstr instead!')
            if nargin == 2
                param = varargin{1};
                if isfield(param,'nExt')
                    param.nExt = param.nExt';
                end
            end
            obj = genmagstr(obj, param);
        end
        function [obj, stat] = sw_anneal(obj, varargin)
            warning('sw:sw_anneal:Obsolete','sw_anneal is obsolete, use anneal instead!')
            [obj, stat] = anneal(obj, varargin{:});
        end
        function E = sw_e(obj, varargin)
            warning('sw:sw_e:Obsolete','sw_e is obsolete, use energy instead!')
            E = energy(obj, varargin{:});
        end
        function F2 = sw_fsf(obj, varargin)
            warning('sw:sw_fsf:Obsolete','sw_fsf is obsolete, use structfact instead!')
            F2 = structfact(obj, varargin{:});
        end
        function spectraTri = sw_spinwave(hTri,hkl,varargin)
            warning('sw:sw_spinwave:Obsolete','sw_spinwave is obsolete, use spinwave instead!')
            spectraTri = spinwave(hTri,hkl,varargin{:});
        end
        function [obj, x, e, exitflag, output] = sw_optmagstr(obj, varargin)
            warning('sw:sw_optmagstr:Obsolete','sw_optmagstr is obsolete, use optmagstr instead!')
            [obj, x, e, exitflag, output] = optmagstr(obj, varargin{:});
        end
    end % classdef
    
end