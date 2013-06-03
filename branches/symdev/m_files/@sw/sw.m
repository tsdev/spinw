classdef sw < class_handlelight
    % SW class defines data structure and methods to calculate spin wave
    % dispersion in magnetic crystals.
    %
    % SW(obj) constructs new sw class object. If obj is sw class, it only
    % checks its data integrity. If obj is struct type, it creates new sw
    % object and checks data integrity.
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
    % Documentation can be found here:
    % <a href='https://wiki.helmholtz-berlin.de/spinw'>https://wiki.helmholtz-berlin.de/spinw</a>
    % Forum for questions:
    % <a href='https://groups.google.com/forum/#!forum/spinwforum'>https://groups.google.com/forum/#!forum/spinwforum</a>
    % Lates version and bug reports/feature requests:
    % <a href='http://code.google.com/p/spinw/'>http://code.google.com/p/spinw/</a>
    %
    %
    
    properties
        lattice     % lattice parameters, fields: angle, lat_const, sym
        unit_cell   % atoms in the unit cell, fields: r, S, label, color
        twin        % twins, fields: vol, rotc 
        matrix      % definition of 3x3 matrices, fields: mat, color, label
        single_ion  % single ion terms in the Hamiltonian, fields: aniso, field
        coupling    % magnetic interactions, fields: dl, atom1, atom2, mat_idx, idx
        mag_str     % magnetic structure, fields: N_ext, k, S, n
        unit        % units of energy, magnetic field and temperature
    end
    
    properties (Access = private)
        matomstore = [];
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
                fprintf(1,'Crystal structure is imported from the given .cif file.\n');
                objS = initfield(struct);
                
                cif0 = cif(firstArg);
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
                    %objS.unit_cell.S = zeros(1,nAtom);
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
