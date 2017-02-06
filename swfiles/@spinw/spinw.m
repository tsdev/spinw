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
    % obj = SPINW(file_path)
    %
    % construct new spinw class object, where file_path contains a string
    % of a .cif/.fst file defining an input crystal and/or magnetic
    % structure.
    %
    % obj = SPINW(figure_handle)
    %
    % copy the SpinW object stored in a previous structural plot.
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
    % spinw Methods:
    % Lattice operations:
    %   genlattice      - defines lattice parameters, angles and symmetry
    %   basisvector     - returns the basis vectors in a matrix
    %   rl              - returns the reciprocal lattice vectors in a matrix
    %   nosym           - reduces the lattice symmetry to P0
    %   newcell         - generates an arbitrary new cell and reduces symmetry to P0
    %   addatom         - adds a new atoms to the lattice
    %   unitcell        - selects a subset of the symmetry inequivalent atoms 
    %   abc             - returns a vector [a b c alpha beta gamma]
    %   atom            - returns all symmetry generated atom
    %   matom           - returns all symmetry generated magnetic atoms
    %   natom           - returns the number of symmetry inequivalent atoms
    %   formula         - returns basic unit cell information
    %   disp            - return basic lattice information
    %   symmetry        - returns whether symmetry is on (>P0)
    %   
    % Plotting:
    %   plot            - plots crystal, magnetic structure and spin Hamiltonian
    %
    % Crystallographic twin operations:
    %   addtwin         - adds crystallographic twin
    %   twinq           - determines reciprocal lattice positions in twins
    %   notwin          - removes the twins
    %   ntwin           - number of twins
    %
    % Magnetic structure operations:
    %   genmagstr       - generate different types of magnetic structures
    %   magstr          - returns the magnetic structure in a rotating frame representation
    %   magtable        - returns a table of the magnetic structure
    %   nmagext         - returns the number of magnetic moments in the supercell
    %   optmagstr       - optimizes magnetic structure with constraints
    %   optmagk         - optimizes the magnetic propagation vector
    %   optmagsteep     - optimizes the moment directions for fixed-k
    %   anneal          - simulated annealing of magnetic structures
    %   annealloop      - simulated annealing and scanning a parameter in the Hamiltonian
    %   structfact      - calculates magnetic structure factor (EXPERIMENTAL)
    %   
    % Matrix operations:
    %   addmatrix        - add a new matrix to the internal data
    %   getmatrix       - determines the symmetry allowed elements of a matrix
    %   setmatrix       - generates matrix values based on symmetry allowed elements
    %   nmat            - number of matrices
    %   
    % Spin Hamiltonian generations:
    %   quickham        - quick generation of Heisenberg Hamiltonian
    %   gencoupling     - generates the list of bonds
    %   addcoupling     - assigns a matrix to a bond
    %   couplingtable   - lists information on the bonds
    %   addaniso        - assigns a matrix to a magnetic ion as anisotropy
    %   addg            - assigns a matrix to a magnetic ion as g-tensor
    %   field           - stores the magnetic field
    %   nbond           - number of bonds
    %   temperature     - temperature for thermal population calculation
    %   intmatrix       - returns the spin Hamiltonian after symmetry applied
    %   
    % Calculators:
    %   spinwave        - linear spin wave solver
    %   powspec         - powder spectrum calculator
    %   energy          - ground state energy
    %   moment          - moment reduction due to quantum fluctuations
    %   spinwavesym     - symbolic spin wave solver
    %   symbolic        - returns whether symbolic mode is on
    %   meanfield       - mean field calculation of q-dependent susceptibility (EXPERIMENTAL)
    %   fourier         - fourier transformation of the Hamiltonian (EXPERIMENTAL)
    %   fouriersym      - symbolic fourier transformation (EXPERIMENTAL)
    %
    % Fitting spin wave spectrum:
    %   fitspec         - fits spin wave energies using global optimization
    %   matparser       - assigns new matrix values based on selectors
    %   horace          - outputs spectrum for Horace
    %   
    % Miscellaneous:
    %   copy            - hard copy of SpinW object
    %   export          - exports magnetic structure into file
    %   fileid          - controls text output of commands
    %   table           - formatted output of internal data
    %   validate        - validates SpinW object
    %   version         - version of SpinW used to create the object
    %   struct          - convert SpinW ojbect to struct
    %   spinw           - constructor
    %
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
        %   lat_const   lattice constants in a 1x3 vector in Angstrom units
        %   angle       (alpha,beta,gamma) angles in 1x3 vector in radian
        %   sym         crystal space group, line number in symmetry.dat file
        %   origin      Origin of the cell in l.u.
        %   label       Label of the space group.
        %
        % See also SPINW.GENLATTICE, SPINW.ABC, SPINW.BASISVECTOR, SPINW.NOSYM.
        lattice
        % Stores the atoms in the crystallographic unit cell.
        % Sub fields are:
        %   r       positions of the atoms in the unit cell, in a
        %           3 x nAtom matrix, in lattice units
        %   S       spin quantum number of the atoms, in a 1 x nAtom
        %           vector, non-magnetic atoms have S=0
        %   label   label of the atom, strings in a 1 x nAtom cell
        %   color   color of the atom in 3 x nAtom matrix, where every
        %           column is an 0-255 RGB color
        %   ox      oxidation number of the atom, in a 1 x nAtom matrix
        %   occ     site occupancy in a 1 x nAtom matrix
        %   b       scattering length of the site for neutron and x-ray
        %           stored in a 2 x nAtom matrix, first row is neutron,
        %           second row is for x-ray
        %   ff      form factor of the site stored in a 2 x 9 x nAtom
        %           matrix, first row is the magnetic form factor for
        %           neutrons, the second row is the charge form factor
        %           for x-ray cross section
        %   Z       atomic number
        %   A       atomic mass (N+Z) for isotopes and -1 for natural
        %           distribution of isotopes
        %   biso    Isotropic displacement factors in units of Angstrom^2.
        %           Definition is the same as in FullProf, defining the
        %           Debye-Waller factor as:
        %               Wd = 1/8*biso/d^2
        %           including in the structure factor as exp(-2Wd)
        %
        % See also SPINW.ADDATOM, SPINW.ATOM, SPINW.MATOM, SPINW.NEWCELL, SPINW.PLOT.
        unit_cell
        % Stores the crystallographic twin parameters.
        % Sub fields are:
        %   rotc    rotation matrices in the xyz coordinate system for
        %           every twin, stored in a 3 x 3 x nTwin matrix
        %   vol     volume ratio of the different twins, stored in a
        %           1 x nTwin vector
        %
        % See also SPINW.ADDTWIN, SPINW.TWINQ, SPINW.UNIT_CELL.
        twin
        % Stores 3x3 matrices for using them in the Hailtonian.
        % Sub fields are:
        %   mat     stores the actual values of 3x3 matrices, in a
        %           3 x 3 x nMatrix matrix, defult unit is meV
        %   color   color assigned for every matrix, stored in a
        %           3 x nMatrix matrix, with 0-255 RGB columns
        %   label   label for every matrix, stored as string in a
        %           1 x nMatrix cell
        %
        % See also SPINW.ADDMATRIX, SPINW.NTWIN.
        matrix
        % Stores single ion terms of the Hamiltonian.
        % Sub fields are:
        %   aniso   vector contains 1 x nMagAtom integers, each integer
        %           assignes one of the nMatrix from the .matrix field
        %           to a magnetic atom in the spinw.matom list as a single
        %           ion anisotropy (zeros for no anisotropy)
        %   g       vector contains 1 x nMagAtom integers, each integer
        %           assignes one of the nMatrix from the .matrix field
        %           to a magnetic atom in the spinw.matom list as a
        %           g-tensor
        %   field   external magnetic field stored in a 1x3 vector,
        %           default unit is Tesla
        %   T       temperature, scalar, default unit is Kelvin
        %
        % See also SPINW.ADDANISO, SPINW.ADDG, SPINW.GETMATRIX, SPINW.SETMATRIX, SPINW.INTMATRIX.
        single_ion
        % Stores the list of spin-spin couplings.
        % Sub fields are:
        %   dl      distance between the unit cells of two interacting
        %           spins, stored in a 3 x nCoupling matrix
        %   atom1   first magnetic atom, pointing to the list of
        %           magnetic atoms in spinw.matom list, stored in a
        %           1 x nCoupling vector
        %   atom2   second magnetic atom, stored in a  1 x nCoupling
        %           vector
        %   mat_idx stores pointers to matrices for every coupling in a
        %           3 x nCoupling matrix, maximum three matrices per
        %           coupling (zeros for no coupling)
        %   idx     increasing indices for the symmetry equivalent
        %           couplings, starting with 1,2,3...
        %   type    Type of coupling corresponding to mat_idx matrices.
        %           Default is 0 for quadratic exchange. type = 1 for
        %           biquadratic exchange.
        %   sym     If true the bond symmetry operators will be applied
        %           on the assigned matrix.
        %   rdip    Maximum distance until the dipolar interaction is
        %           calculated. Zero vlue means no dipolar interactions
        %           are considered.
        %   nsym    The largest bond 'idx' value that is generated
        %           using the space group operators. Typically very
        %           long bonds for dipolar interactions won't be
        %           calculated using space group symmetry.
        %
        % See also SPINW.GENCOUPLING, SPINW.ADDCOUPLING, SPINW.FIELD.
        coupling
        % Stores the magnetic structure.
        % Sub fields are:
        %   S       stores the moment direction for every spin in the
        %           crystallographic or magnetic supercell in a
        %           3 x nMagExt matrix, where nMagExt = nMagAtom*prod(N_ext)
        %   k       magnetic ordering wave vector in a 3x1 vector
        %   n       normal vector to the rotation of the moments in
        %           case of non-zero ordering wave vector, dimensions are
        %           3x1
        %   N_ext   Size of the magnetic supercell, default is [1 1 1]
        %           if the magnetic cell is identical to the
        %           crystallographic cell, the 1x3 vector extends the
        %           cell along the a, b and c axis
        %
        % See also SPINW.GENMAGSTR, SPINW.OPTMAGSTR, SPINW.ANNEAL, SPINW.MOMENT, SPINW.NMAGEXT, SPINW.STRUCTFACT.
        mag_str
        % Stores the physical units in the Hamiltonian.
        % Defaults are meV, Tesla Angstrom and Kelvin.
        % Sub fields are:
        %   kB      Boltzmann constant, default is 0.0862 [meV/K]
        %   muB     Bohr magneton, default is 0.0579 [meV/T]
        %   mu0     vacuum permeability, 201.335431 [T^2*Angstrom^3/meV]
        unit
        % Stores the cache, it should be only used to check consistency of the code.
        % The stored values should not be changed by the user in any case!
        % Sub fields are:
        %   matom   data on the magnetic unit cell
        %   symop   data on the generated symmetry operators per bond
        cache = struct('matom',[],'symop',[]);
        
    end
    
    properties (Access = private)
        % stores the property change listener handles
        propl = event.proplistener.empty;
        % stores whether the couplings are generated under symmetry constraints
        sym   = false;
        % stores whether the calculation are done symbolically
        symb  = false;
        % stores the file ID of the text output, default is the Command Window (see swpref)
        fid   = 1;
        % use the version property as contant, this will be executed only
        % once
        ver   = sw_version;
    end
    
    methods
        function obj = spinw(varargin)
            % SPINW constructor
            %
            
            % update fid value
            obj.fid = swpref.getpref('fid',[]);

            if nargin==0
                objS = initfield(struct);
                fNames = fieldnames(objS);
                for ii = 1:length(fNames)
                   obj.(fNames{ii}) = objS.(fNames{ii});
                end
                return
            end
            
            firstArg = varargin{1};
            
            if ishandle(firstArg)
                % get spinw object from graphics handle
                switch get(firstArg,'Tag')
                    case 'sw_crystal'
                        figDat = getappdata(firstArg);
                        obj = copy(figDat.obj);
                    case 'sw_spectra'
                        figDat = getappdata(firstArg);
                        obj    = copy(figDat.spectra.obj);
                end
                return

            end
            
            if isa(firstArg, 'spinw')
                %  it is used when objects are passed as arguments.
                obj = copy(firstArg);
                return
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
                % import data from file (cif/fst are supported)
                
                objS = initfield(struct);
                fNames = fieldnames(objS);
                for ii = 1:length(fNames)
                    obj.(fNames{ii}) = objS.(fNames{ii});
                end
                
                obj = sw_import(firstArg,false,obj);
                
            end
            
        end % .spinw
        
        
        function nMagExt = nmagext(obj)
            % gives the number of magnetic atoms in the magnetic supercell
            %
            % nMagExt = NMAGEXT(obj)
            %
            
            nMagExt = size(obj.mag_str.F,2);
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
    
    methods(Hidden=true,Static=true)
            function obj = loadobj(obj)
                % restore property listeners
                % add new listeners to the new object
                if ~isempty(obj.cache.matom)
                    % add listener to lattice and unit_cell fields
                    obj.addlistenermulti(1);
                end
                if ~isempty(obj.cache.symop)
                    % add listener to lattice, unit_cell and coupling fields
                    obj.addlistenermulti(2);
                end
            end
    end
        
    methods(Hidden=true)
        function clearcache(obj, chgField)
            % clears the cache
            %
            % CLEARCACHE(obj)
            %
            
            % listening to changes of the spinw object to clear cache is
            % necessary
            
            if nargin<2
                % delete the existing listener handles
                delete(obj.propl(ishandle(obj.propl)));
                % remove cache
                obj.cache.matom = [];
                obj.cache.symop = [];
                return
            end
            
            switch chgField
                case 1
                    % magnetic atoms: delete the stored magnetic atom positions
                    obj.cache.matom = [];
                    % remove the listeners
                    delete(obj.propl(1:2));
                case 2
                    % bond symmetry operators: delete the stored operators
                    obj.cache.symop = [];
                    % remove the listeners
                    delete(obj.propl(3:5));
            end
        end
        
        function addlistenermulti(obj, chgField)
            % create the corresponding listeners to each cache subfield
            
            switch chgField
                case 1
                    % add listener to lattice and unit_cell fields
                    obj.propl(1) = addlistener(obj,'lattice',  'PostSet',@(evnt,src)obj.clearcache(1));
                    obj.propl(2) = addlistener(obj,'unit_cell','PostSet',@(evnt,src)obj.clearcache(1));
                case 2
                    % add listener to lattice, unit_cell and coupling fields
                    obj.propl(3) = addlistener(obj,'lattice',  'PostSet',@(evnt,src)obj.clearcache(2));
                    obj.propl(4) = addlistener(obj,'unit_cell','PostSet',@(evnt,src)obj.clearcache(2));
                    obj.propl(5) = addlistener(obj,'coupling', 'PostSet',@(evnt,src)obj.clearcache(2));
            end
        end
        
        function obj = saveobj(obj)
            % remove property change listeners
            delete(obj.propl);
            % empty pointers
            obj.propl = event.listener.empty;
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