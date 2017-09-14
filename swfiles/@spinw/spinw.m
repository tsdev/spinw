classdef spinw < handle & matlab.mixin.SetGet
    % class to store and solve magnetic Hamiltonians
    %
    % * * *
    % `obj = spinw()`
    % * * *
    %
    % constructs a new spinw class object, with default parameters.
    %
    % * * *
    % `obj = spinw(obj)`
    % * * *
    %
    % constructs new spinw class object. If `obj` is spinw class, it only
    % checks the integrity of its internal data structure. If `obj` is
    % struct type, it creates new spinw object and checks data integrity.
    %
    % * * *
    % `obj = spinw(source)`
    % * * *
    %
    % construct new spinw class object, where `source` is either a file
    % path pointing to a local cif or fst file or a link to an online file.
    %
    % * * *
    % `obj = spinw(figure_handle)`
    % * * *
    %
    % copy the spinw object stored in a previous structural3D plot figure.
    %
    % The data structure behind the spinw object can be accessed by using
    % `struct(obj)`. All fields of the struct type data behind the spinw
    % object are accessible through the main field names of the `obj`
    % object. For example the lattice parameters can be extracted using:
    % ```matlab
    %   abc = obj.unit_cell.lat_const
    % ```
    %
    % spinw is a handle class, that means that only the handle of the
    % object is copied in a `swobj1 = swobj2` command. To create a copy
    % (clone) of an spinw object use
    % ```matlab
    %    swobj1 = swobj2.copy
    % ```
    %
    % ### Methods
    %
    % #### Lattice operations
    %
    %   spinw.genlattice
    %   spinw.basisvector
    %   spinw.rl
    %   spinw.nosym
    %   spinw.newcell
    %   spinw.addatom
    %   spinw.unitcell
    %   spinw.abc
    %   spinw.atom
    %   spinw.matom
    %   spinw.natom
    %   spinw.formula
    %   spinw.disp
    %   spinw.symmetry
    %   
    % #### Plotting
    %
    %   spinw.plot
    %
    % #### Crystallographic twin operations
    %
    %   spinw.addtwin
    %   spinw.twinq
    %   spinw.notwin
    %   spinw.ntwin
    %
    % #### Magnetic structure operations
    %
    %   spinw.genmagstr
    %   spinw.magstr
    %   spinw.magtable
    %   spinw.nmagext
    %   spinw.optmagstr
    %   spinw.optmagk
    %   spinw.optmagsteep
    %   spinw.anneal
    %   spinw.annealloop
    %   spinw.structfact
    %   
    % #### Matrix operations
    %
    %   spinw.addmatrix
    %   spinw.getmatrix
    %   spinw.setmatrix
    %   spinw.nmat
    %   
    % #### Spin Hamiltonian generations
    %
    %   spinw.quickham
    %   spinw.gencoupling
    %   spinw.addcoupling
    %   spinw.couplingtable
    %   spinw.addaniso
    %   spinw.addg
    %   spinw.field
    %   spinw.nbond
    %   spinw.temperature
    %   spinw.intmatrix
    %   spinw.symop
    %   spinw.setunit
    %   
    % #### Calculators
    %
    %   spinw.spinwave
    %   spinw.powspec
    %   spinw.energy
    %   spinw.moment
    %   spinw.spinwavesym
    %   spinw.symbolic
    %   spinw.meanfield
    %   spinw.fourier
    %   spinw.fouriersym
    %
    % #### Fitting spin wave spectrum
    %
    %   spinw.fitspec
    %   spinw.matparser
    %   spinw.horace
    %   
    % #### Miscellaneous
    %
    %   spinw.copy
    %   spinw.export
    %   spinw.fileid
    %   spinw.table
    %   spinw.validate
    %   spinw.version
    %   spinw.struct
    %   spinw.clearcache
    %   spinw.spinw
    %
    % ### Notes
    %
    % Tutorials and documentation can be found at [psi.ch/spinw](https://psi.ch/spinw)
    %
    % Forum for questions on [Google Groups](https://groups.google.com/forum/#!forum/spinwforum)
    %
    % Lates version and bug reports/feature requests can be submitted on [GitHub](https://github.com/tsdev/spinw)
    %
    % ### See also
    %
    % [spinw.copy], [spinw.struct], [Comparing handle and value classes](https://www.google.ch/url?sa=t&rct=j&q=&esrc=s&source=web&cd=3&ved=0ahUKEwjCvbbctqTWAhVBblAKHQxnAnIQFggyMAI&url=https%3A%2F%2Fwww.mathworks.com%2Fhelp%2Fmatlab%2Fmatlab_oop%2Fcomparing-handle-and-value-classes.html&usg=AFQjCNFoN4qQdn6rPXKWkQ7aoog9G-nHgA)
    %

    properties (SetObservable)
        % stores the unit cell parameters
        %
        % ### Sub fields
        %
        % `lat_const`
        % : Lattice constants in a $[1\times 3]$ vector in units defined in
        %   [spinw.unit] (default value is \\Angstrom).
        %
        % `angle`
        % : `[\\alpha,\\beta,\\gamma]` angles in a $[1\times 3]$ vector in
        %   radian units.
        %
        % `sym`
        % : Symmetry operators stored in matrix with dimensions of
        %   $[3\times 4 \times n_{op}]$.
        %
        % `origin`
        % : Origin of the cell in lattice units.
        % 
        % `label`
        % : Label of the space group.
        %
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
        % stores single ion terms of the Hamiltonian
        %
        % ### Sub fields
        %
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
        % stores the list of bonds
        %
        % ### Sub fields
        %
        % `dl`
        % : Distance between the unit cells of two interacting
        %   spins, stored in a $[3\times n_{coupling}]$ matrix.
        %
        % `atom1`
        % : First magnetic atom, pointing to the list of
        %   magnetic atoms in the list generated by `spinw.matom`, stored in a
        %   row vector with $n_{coupling}$ element.
        %
        % `atom2`
        % : Second magnetic atom, same as `atom1`.
        %
        % `mat_idx`
        % : Stores pointers to matrices for every coupling in a
        %   $[3\times n_{coupling}]$ matrix, maximum three matrix per
        %   coupling (zeros for no coupling) is allowed.
        %
        % `idx`
        % : Neighbor index, increasing indices for the equivalent
        %   couplings, starting with 1,2,... which means first and second
        %   neighbor bonds, respectively.
        %
        % `type`
        % : Type of coupling corresponding to `mat_idx` matrices.
        %   Default is 0 for quadratic exchange, `type = 1` for
        %   biquadratic exchange.
        % 
        % `sym`
        % : If `true` the bond symmetry operators will be applied
        %   on the assigned matrix.
        % 
        % `rdip`
        % : Maximum distance until the dipolar interaction is
        %   calculated. Zero value means no dipolar interactions
        %   are considered.
        % 
        % `nsym`
        % : The largest bond `idx` value that is generated
        %   using the space group operators. Typically very long bonds for
        %   dipolar interactions won't be calculated using space group
        %   symmetry.
        %
        coupling
        % stores the magnetic structure
        %
        % ### Sub fields
        %
        % `F`
        % : Complex magnetization (strictly speaking complex
        %   spin expectation value) for every spin in the magnetic
        %   cell, represented by a matrix with dimensions of $[3\times
        %   n_{magext}\times n_k]$,
        %   where `nMagExt = nMagAtom*prod(N_ext)` and $n_k$ is the number
        %   of the magnetic propagation vectors.
        %
        % `k`
        % : Magnetic propagation vectors stored in a matrix with dimensions
        %   of $[3\times n_k]$.
        %
        % `N_ext`
        % : Size of the magnetic supercell in lattice units, default value
        %   is `[1 1 1]` emaning that the magnetic cell is identical to the
        %   crystallographic cell. The $[1\times 3]$ vector extends the cell
        %   along the $a$, $b$ and $c$ axes.
        %
        mag_str
        % stores the physical units for the Hamiltonian
        %
        % Default values are meV, T, \\Angstrom and K for energy, magnetic
        % field, length and temperature, respectively.
        %
        % ### Sub fields
        %
        % `kB`
        % : Boltzmann constant, default value is 0.0862 meV/K.
        %
        % `muB`
        % : Bohr magneton, default values is 0.0579 meV/T.
        %
        % `mu0`
        % : Vacuum permeability, default value is 201.335431 T$^2$\\Angstrom$^3$/meV.
        %
        % `label`
        % : Labels for distance, energy, magnetic field and temperature
        % stored in a cell with dimensions $[1\times 4]$.
        %
        % `nformula`
        % : Number of formula units in the unit cell.
        %
        % `qmat`
        % : Transformation matrix that converts the given $Q$ values with
        % dimensions of $[3\times 3]$.
        %
        unit
        % stores temporary values
        %
        % This property should be only used to check consistency of the code.
        %
        % {% include warning.html content="Changing these values is strongly 
        % discouraged as it can break the code!" %}
        % 
        % ### Sub fields
        %
        % `matom`
        % : Generated data of the magnetic unit cell.
        %
        % `symop`
        % : Generated symmetry operators for each bond.
        %
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
    
    methods (Static)
        % static methods
        validate(varargin)
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
                
                spinw.validate(objS);
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
            obj.propl = event.proplistener.empty;
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