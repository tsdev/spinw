function varargout = formula(obj, print)
% returns chemical formula, mass, volume, etc.
%
% formula = FORMULA(obj,{print})
%
% Options:
%
% obj       sw class object.
% print     If true, the results are printed onto the Matlab Command
%           Window. Default is false. Optional
%
% Output:
%
% formula struct variable with the following fields:
% m         Mass of the unit cell in g/mol unit.
% V         Volume of the unit cell in Angstrom^3 unit.
% rho       Density in g/cm^3 unit.
% chemlabel List of the different elements.
% chemnum   Number of the listed element names
% chemform  Chemical formula string: series of 'ChemLabel_ChemNum '.
%
% Example:
%
% cryst = sw('test.cif')
% cryst.formula(true);
%
% The formula of the crystal stored in the test.cif file will be printed
% onto the Command Window.
%

if nargin == 1
    print = false;
end

atom = obj.atom;

m = zeros(1,numel(atom.idx));
aLabel = cell(1,numel(atom.idx));

for ii = 1:numel(atom.idx)
    [m(ii), aLabel{ii}] = sw_atomdata(obj.unit_cell.label{atom.idx(ii)},'mass');
end

if numel(atom.idx) == 0
    m = 0;
end

% aVogadro number (1/mol)
nA = 6.02214129e23;

% mass in gramms
formula.m = sum(m);
% volume in Angstrom^3
formula.V = prod(obj.lattice.lat_const);
% density in g/cm^3
formula.rho = formula.m/formula.V/nA*1e24;

% find fomula
diffLabel = unique(aLabel);

numAtom = cell(1,2*numel(diffLabel));

for ii = 1:numel(diffLabel)
    temp = sum(strcmp(aLabel,diffLabel{ii}));
    numAtom{2*ii-1} = diffLabel{ii};
    numAtom{2*ii} = temp;
end

formula.chemform = sprintf('%s_%d ',numAtom{:});
formula.chemlabel = numAtom(1:2:end);
formula.chemnum = [numAtom{2:2:end}];

if nargout > 0
    varargout{1} = formula;
end

if print
    fprintf('Chemical formula: %s\n',formula.chemform);
    fprintf('Mass:             %8.3f g/mol\n',formula.m);
    fprintf('Volume:           %8.3f Angstrom^3\n',formula.V);
    fprintf('Density:          %8.3f g/cm^3\n',formula.rho);
    
end

end