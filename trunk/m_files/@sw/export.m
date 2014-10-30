function out = export(obj, varargin)
% export data from sw object into different file formats
%
% out = EXPORT(obj, 'Option1', Value1, ...)
%
% Different part of the sw object data can be exported selected by the
% 'format' option. Right now the following formats are supported:
%
% 'pcr'     Creates part of a .pcr file used by FullProf. It exports the
%           atomic positions.
% 'spt'     Creates a Jmol script, that reproduce the same plot as used the
%           built in sw.plot() function. Any additional parameter of the
%           sw.plot() function can be used.
% 'MC'      Exports a proprietary file format for Monte Carlo simulations.
%
%
% Other general options:
%
% path      Path to a file into which the data will be exported, 'out' will
%           be true if the file succesfully saved, otherwise false.
% fid       File identifier that is already opened in Matlab using the
%           fid = fopen() function. 'out' will be the input fid. Don't
%           forget to close the text file afterwards.
%
%
% Format related options:
%
% PCR format:
% perm      Permutation of the xyz atomic positions, default is [1 2 3].
%
% MC format:
% boundary  Boundary conditions of the extended unit cell.
%                 'free'  Free, interactions between extedned unit cells are
%                         omitted.
%                 'per'   Periodic, interactions between extended unit cells
%                         are retained.
%             Default is {'per' 'per' 'per'}.
%
% If neither 'path' nor 'fid' is given, the 'out' will be a cell containing
% strings for each line of the text output.
%
% Example:
%
% cryst = sw('test.cif');
% cryst.export('format','pcr','path','test.pcr');
%
% In this example the crystal structure is imported from the test.cif file,
% and the atomic positions are saved into the test.pcr file for FullProf
% refinement (the pcr file needs additional text to work with FullProf).
%
% Links:
% Jmol Wiki: http://wiki.jmol.org/index.php/Main_Page
% FullProf:  https://www.ill.eu/sites/fullprof
%

inpForm.fname  = {'format' 'path' 'fid' 'perm'  'boundary'          };
inpForm.defval = {''       ''      []   [1 2 3] {'per' 'per' 'per'} };
inpForm.size   = {[1 -1]   [1 -2] [1 1] [1 3]   [1 3]               };
inpForm.soft   = {true     true    true false   false               };

if nargin == 1
    varargin{1}.showWarn = false;
elseif nargin>1
    varargin{end+1} = 'showWarn';
    varargin{end+1} = false;
end

param = sw_readparam(inpForm, varargin{:});

% produce the requested output

if isempty(param.path) && isempty(param.fid)
    % dialog to get a filename
    [fName, fDir] = uiputfile({'*.pcr','FullProf file (*.pcr)';'*.spt','Jmol script (*.spt)';'*.*' 'All Files (*.*)'}, 'Select an output filename');
    param.path = [fDir fName];
    if isempty(param.format)
        [~,~,fExt] = fileparts(param.path);
        param.format = fExt(2:end);
    end
end


switch param.format
    case 'pcr'
        % create .pcr text file
        outStr = createpcr(obj, param.perm);
    case 'MC'
        outStr = createmc(obj, param.boundary);
    case 'spt'
        % create Jmol script file
        if nargin == 2
            varargin{1}.showWarn = false;
            varargin{1}.format = 'jmol';
        else
            varargin{end+1} = 'showWarn';
            varargin{end+1} = false;
            varargin{end+1} = 'format';
            varargin{end+1} = 'jmol';
            
        end
        
        outStr = plot(obj, varargin{:});
        
    case ''
        warning('sw:export:NoInput','No ''format'' option was given, no output is produced!');
        out = [];
        return
    otherwise
        error('sw:export:WrongInput','''format'' has to be one of the strings given in the help!');
end

% write into fid file
if ~isempty(param.fid)
    fprintf(param.fid,outStr);
    out = param.fid;
    return
end

% save into path file
if ~isempty(param.path)
    try
        fid = fopen(param.path,'w');
        fprintf(fid,outStr);
        fclose(fid);
        out = true;
    catch
        % file couldn't be saved
        out = false;
    end
    return
end

% provide string output
out = outStr;

end

function out = createpcr(obj, perm)
% CREATEPCR(obj, perm) creates the structural part of a pcr file
% from a .cif file.
%
% This function will create the atomic positions from an sw object in the
% input format for FullProf Rietveld refinement software.
%
% perm  Permutation of the (x,y,z) coordinates.
%

% generate all atoms in the unit cell to count site multiplicities
atoms = obj.atom;
mult = accumarray(atoms.idx',ones(numel(atoms.idx),1));
mult = mult/max(mult);

% output string
out = sprintf('!Atom   Typ       X        Y        Z     Biso       Occ     In Fin N_t Spc /Codes\n');

nAtom = size(obj.unit_cell.r,2);

uc = obj.unit_cell;

% split labels into [aname, alabel]
% aname: name of atom (e.g. 'Cr')
% alabel: label if given (eg' 'MCr3'), otherwise the same as the name of
% the atom
uc.aname = cell(1,nAtom);
uc.alabel = cell(1,nAtom);

for ii = 1:nAtom
    lTemp = strword(uc.label{ii},[1 2],true);
    uc.alabel{ii} = lTemp{1};
    uc.aname{ii} = lTemp{2};
end

% find unique labels for atoms
for ii = 1:nAtom
    uc.ulabel(ii) = ~(sum(strcmp(uc.alabel,uc.alabel{ii}))>1);
end

% sort atoms according to the
idx = 0;
for ii = 1:nAtom
    if ~uc.ulabel(ii)
        % not unique atom labels put extra number
        strT = sprintf('%s%d',uc.alabel{ii},idx);
        idx = idx + 1;
    else
        % no extra numbering
        strT = sprintf('%s',uc.alabel{ii});
    end
    % pad the string to 6 characters with whitespace
    if numel(strT)<6
        strT((end+1):6) = ' ';
    end
    strT = [strT sprintf(' %s',uc.aname{ii})]; %#ok<*AGROW>
    % pad the string to 14 characters with whitespace
    if numel(strT)<13
        strT((end+1):13) = ' ';
    end
    
    out = [out strT sprintf('%9.5f%9.5f%9.5f%9.5f%9.5f%4d%4d%4d%4d\n',uc.r(perm,ii)',0,mult(ii),[0 0 0 0])];
    out = [out sprintf('                  0.00     0.00     0.00     0.00      0.00\n')];
end

end


function outStr = createmc(obj, boundary)

% block1: boundary conditions
block1 = zeros(1,3);
for ii = 1:3
    if strcmp('per',boundary{ii})
        block1(ii) = 1;
    end
end

% block2
block2 = reshape(permute(obj.matrix.mat,[3 2 1]),[],9)';

[SS, SI, RR] = obj.intmatrix('zeroC',false,'plotmode',true);
% RR in lattice units
r = bsxfun(@times,double(obj.mag_str.N_ext'),RR)';
% atom index
idx = (1:size(r,1))';
% number of cells
nCell = prod(double(obj.mag_str.N_ext));
% spin of each atom
spin = repmat(obj.matom.S,[1 nCell])';
% block3
block3 = [idx r spin]';

% block4: anisotropy matrices
block4 = reshape(permute(SI.aniso,[3 2 1]),[],9)';

% remove coupling for free boundary conditions
for ii = 1:3
    if strcmp('free',boundary{ii})
        SS.all(:,SS.all(ii,:)~=0) = [];
    end
end

% Since k_m=(0,0,0) the spins that are coupled to themself contribute with
% a constant self-energy, removing this doesn't change thermodynamical
% behaviour just shifts the zero energy.
SS.all(:, SS.all(4,:)==SS.all(5,:)) = [];

% block5: coupling table
block5 = SS.all([end-1 4 5],:);


outStr = sprintf('# boundary conditions (free = 0, periodic = 1)\n');
outStr = [outStr sprintf('%3d %3d %3d\n',block1)];
outStr = [outStr sprintf('# exchange matrices Jxx, Jxy, Jxz, Jyx, ... [9 double per line]\n')];
outStr = [outStr sprintf([repmat('%7.5f ',[1 9]) '\n'],block2)];
outStr = [outStr sprintf('# atom_idx     r_x   r_y   r_z spin\n')];
outStr = [outStr sprintf('%10d %7.3f %5.3f %5.3f %4d\n',block3)];
outStr = [outStr sprintf('# anisotropy matrices Axx, Axy, Axz, Ayx, ... [9 double per line]\n')];
outStr = [outStr sprintf([repmat('%7.5f ',[1 9]) '\n'],block4)];
outStr = [outStr sprintf('# coupling table J_idx atom_idx1 atom_idx2\n')];
outStr = [outStr sprintf('%22d %9d %9d\n',block5)];

end