classdef sw < class_handlelight
    % SW class defines data structure and methods to calculate spin wave
    % dispersion in magnetic crystals.
    %
    % SW() constructs a new sw class object, with default parameters.
    %
    % SW(obj) constructs new sw class object. If obj is sw class, it only
    % checks its data integrity. If obj is struct type, it creates new sw
    % object and checks data integrity.
    %
    % SW(cif_path) construct new sw class object, where cif_path contains a
    % string of a .cif file path defining a crystals structure.
    %
    % The data structure behind the sw object can be accessed by
    % struct(sw). All fields of the struct type data behind the sw object
    % are accessible through the main field names of the sw object. For
    % example the lattice parameters:
    %   abc = sw.unit_cell.lat_const;
    %
    % sw is a handle class, it means that only the handle of the object
    % is copied in a swobj1 = swobj2 command. To create a copy (clone) of
    % an sw object use:
    %    swobj1 = swobj2.copy;
    % See also:
    % <a href='/Applications/MATLAB_R2012b.app/help/matlab/matlab_oop/comparing-handle-and-value-classes.html'>Comparing handle and value classes</a>
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
        %               defalult unit is Tesla
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
        % See also SW.GENMAGSTR, SW.OTPMAGSTR, SW.ANNEAL, SW.MOMENT, SW.NMAGEXT, SW.STRUCTFACT.
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
        issym  = true;  % stores whether the couplings are generated under symmetr constraints
        Elabel = 'meV';
        Qlabel = 'Angstrom^{-1}';
        Rlabel = 'Angstrom';
        Blabel = 'T';
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
                fprintf(1,'Crystal structure is imported from %s.\n',firstArg);
                abc0 = [cif0.cell_length_a cif0.cell_length_b cif0.cell_length_c];
                ang0 = [cif0.cell_angle_alpha cif0.cell_angle_beta cif0.cell_angle_gamma];
                sym0 = cif0.('symmetry_space_group_name_H-M');
                %symi0 = cif0.symmetry_Int_Tables_number;
                xyz0 = cif0.symmetry_equiv_pos_as_xyz; xyz0 = sprintf('%s; ',xyz0{:}); xyz0 = xyz0(1:end-2);
                name0 = cif0.atom_site_type_symbol';
                r0 = [cif0.atom_site_fract_x cif0.atom_site_fract_y cif0.atom_site_fract_z]';
                
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
                        col0 = sw_atomdata(name0{ii}(name0{ii}>57),'color')';
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
            % clones sw object
            %
            
            objS = struct(obj);
            objC = sw(objS);
            
            % copy the private properties
            objC.issym  = obj.issym;
            objC.Elabel = obj.Elabel;
            objC.Qlabel = obj.Qlabel;
            objC.Rlabel = obj.Rlabel;
            objC.Blabel = obj.Blabel;
        end % copy
        function abc = abc(obj)
            % returns [a, b, c, alpha, beta, gamma] vector
            % in Angstrom and degree units
            abc = [obj.lattice.lat_const obj.lattice.angle*180/pi];
        end
        function nMagExt = nmagext(obj)
            % gives the number of magnetic atoms in the extended unit cell
            nMagExt = size(obj.mag_str.S,2);
        end
        function nTwin = ntwin(obj)
            % gives the number of twins
            nTwin = size(obj.twin.vol,2);
        end
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