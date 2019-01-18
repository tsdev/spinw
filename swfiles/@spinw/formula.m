function varargout = formula(obj)
% returns basic physical properties
%
% ### Syntax
%
% `formula = formula(obj)`
%
% ### Description
%
% `result = formula(obj)` returns chemical mass, density, cellvolume etc.
% of `obj`.
%
% ### Examples
%
% The formula of the crystal stored in the
% [https://goo.gl/do6oTh](https://goo.gl/do6oTh) linked file will be
% printed onto the Command Window.
%
% ```
% >>cryst = spinw('https://goo.gl/do6oTh')
% >>cryst.formula>>
% ```
%
% ### Name-Value Pair Arguments
%
% `'obj'`
% : [spinw] object.
%
% ### Output Arguments
%
% `formula` struct variable with the following fields:
% * `m`         Mass of the unit cell in g/mol units.
% * `V`         Calculated volume of the unit cell in length units (defined in [spinw.unit]).
% * `rho`       Density in g/cm$^3$.
% * `chemlabel` List of the different elements.
% * `chemnum`   Number of the listed element names
% * `chemform`  Chemical formula string.
%

atom = obj.atom;

[m, aLabel] = sw_atomdata(obj.unit_cell.label(atom.idx),'mass');

% get more precise mass from isotope data
iso0  = sw_readtable([sw_rootdir 'dat_files' filesep 'isotope.dat']);
iso.m = [iso0(:).mass];
iso.Z = [iso0(:).Z];
iso.A = [iso0(:).A];

% occupancy
occ    = obj.unit_cell.occ(atom.idx);

if numel(atom.idx) == 0
    m = 0;
else
    formulaA = obj.unit_cell.A(atom.idx);
    formulaZ = obj.unit_cell.Z(atom.idx);
    
    [iFound, loc] = ismember([formulaA;formulaZ]',[iso.A;iso.Z]','rows');
    
    % only keep if all atoms are found
    if all(iFound)
        m = iso.m(loc).*occ;
    end
end

% aVogadro number (1/mol)
nA = 6.02214129e23;

% find formula
diffLabel = unique(aLabel);

numAtom = cell(1,2*numel(diffLabel));

for ii = 1:numel(diffLabel)
    %findAtom = strcmp(aLabel,diffLabel{ii});
    temp = sum(strcmp(aLabel,diffLabel{ii}).*occ);
    numAtom{2*ii-1} = diffLabel{ii};
    numAtom{2*ii} = temp;
end

% number of formula units in the unit cell
nForm = double(obj.unit.nformula);
if nForm == 0
    % formula units is not given directly, try to find greates commond divider
    nForm = gcdv([numAtom{2:2:end}]);
    if nForm < 1
        nForm = 1;
    end
end

% number of formula in cell
formula.N = nForm;
% formula mass in gramms
formula.m = sum(m)/nForm;
% cell volume in Angstrom^3
formula.V = det(obj.basisvector);
% density in g/cm^3 convert volume to A^3
formula.rho = formula.m/(formula.V*sw_converter(1,obj.unit.label{1},symbol('a'))^3)/nA*1e24*nForm;

% divide the number of atoms in formula
numAtom(2:2:end) = num2cell([numAtom{2:2:end}]/nForm);

formula.chemform = sprintf('%s%g',numAtom{:});
formula.chemlabel = numAtom(1:2:end);
formula.chemnum = [numAtom{2:2:end}];

if nargout > 0
    varargout{1} = formula;
end

if nargout == 0
    str = [...
        sprintf('     `Chemical formula:`  %s\n',formula.chemform)...
        sprintf('     `Formula mass:`      %8.3f g/mol\n',formula.m)...
        sprintf('     `Formula in cell:`   %8d units\n',formula.N)...
        sprintf('     `Cell volume:`       %8.3f %s\\\\^3\n',formula.V,obj.unit.label{1})...
        sprintf('     `Density:`           %8.3f g/cm\\\\^3\n',formula.rho)];
    fprintf(sw_markdown(str));
end

end