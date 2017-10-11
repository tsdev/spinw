function varargout = sw_atomdata(atomSymb, datType)
% returns information on chemical elements
% 
% ### Syntax
% 
% `data = sw_atomdata(atomsymb)`
% 
% `data = sw_atomdata(atomsymb,datatype)`
%
% ### Description
% 
% `data = sw_atomdata(atomsymb)` returns information on chemical elements
% (RGB color code, mass, long name) in a struct. The element is identified
% by its short name, such as 'O' for oxygen. If the given atom name does
% not exists, the function returns the data for `'Unobtanium'`.
% 
% `data = sw_atomdata(atomsymb,datatype)` returns only the requested type
% of data.
%
% ### Examples
%
% The radius of the hydrogen atom:
% ```
% >>r_H = sw_atomdata('H','radius')>>
% ```
% 
% ### Input Arguments
% 
% `atomSymb`
% : String of the name of the atom, for example 'He' or the atomic
%   number Z. If the string contains whitespace character, the
%   second word will be used to identify the atom.
% 
% `dataType`
% : Type of information requested, following strings are accepted:
%   * `'name'`        Atomic symbol.
%   * `'radius'`      Atomic radius.
%   * `'color'`       Color of the atom from the [CPK color scheme](https://en.wikipedia.org/wiki/CPK_coloring).
%   * `'mass'`        Average mass of the element.
%   * `'longname'`    Name of the element.
%   * `'Z'`           Atomic index.
%   * `'all'`         All atomic data returned in a struct.
% 
% ### See Also
% 
% [sw_mff] \| [sw_cff]
%

if nargin == 1
    datType = 'all';
end

if nargin == 0
    help sw_atomdata
    return
end

% read in the atom definition file
atom  = sw_readtable([sw_rootdir 'dat_files' filesep 'atom.dat']);

% add atomic number for the first 113 atom + unobtanium (original definitions)
cellIdx = num2cell(1:113);
% FANCY
[atom(1:113).Z] = deal(cellIdx{:});

if ischar(atomSymb)
    atomSymb = {atomSymb};
    backChar = true;
else
    backChar = false;
end

if iscell(atomSymb)
    % point to the empty atom data by default
    idx = zeros(1,numel(atomSymb))+numel(atom);
    
    for ii = 1:numel(atomSymb)
        % split multiple words and use the second word if exists
        atomSymb0 = strword(atomSymb{ii},2,true);
        atomSymb0 = atomSymb0{1};
        
        % cut M from the beginning of the atom label
        upStr = isstrprop(atomSymb0,'upper');
        if (numel(atomSymb0)>=2) && all(upStr(1:2))
            atomSymb0 = atomSymb0(2:end);
        end
        % cut end symbols
        atomSymb0 = atomSymb0(isstrprop(atomSymb0,'alpha'));
        % find atom symbol
        idx0 = find(strcmp({atom(:).name},atomSymb0));
        
        atomSymb{ii} = atomSymb0;
        
        if ~isempty(idx0)
            idx(ii) = idx0;
        else
            atomSymb{ii} = 'A';
        end
    end
    
    if backChar
        atomSymb = atomSymb{1};
    end
else
    idx = atomSymb;
    idx(idx<1 | idx>113) = 113;
end

% select atoms
atom = atom(idx);

switch datType
    case 'name'
        data = {atom.name};
    case 'radius'
        data = [atom.R];
    case 'color'
        data = reshape([atom.RGB],3,[]);
    case 'mass'
        data = [atom.mass];
    case 'longname'
        data = {atom.longname};
    case 'Z'
        data = [atom.Z];
    case 'all'
        data = rmfield(atom,'MODE');
    otherwise
        error('sw_atomdata:WrongInput','datType has to be one of the string options, see help sw_atomdata!');
end

varargout{1} = data;

if nargout > 1
    varargout{2} = atomSymb;
end

end