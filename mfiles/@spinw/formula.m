function varargout = formula(obj,fid)
% returns chemical formula, mass, volume, etc.
%
% formula = FORMULA(obj)
%
% Options:
%
% obj       spinw class object.
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
% cryst.formula;
%
% The formula of the crystal stored in the test.cif file will be printed
% onto the Command Window.
%

if nargin < 2
    fid = obj.fid;
end

atom = obj.atom;

[m, aLabel] = sw_atomdata(obj.unit_cell.label(atom.idx),'mass');

% get more precise mass from isotope data
iso0  = sw_readtable([sw_rootdir 'dat_files' filesep 'isotope.dat']);
iso.m = [iso0(:).mass];
iso.Z = [iso0(:).Z];
iso.A = [iso0(:).A];

if numel(atom.idx) == 0
    m = 0;
else
    formulaA = obj.unit_cell.A(atom.idx);
    formulaZ = obj.unit_cell.Z(atom.idx);
    
    [iFound, loc] = ismember([formulaA;formulaZ]',[iso.A;iso.Z]','rows');
    
    % only keep if all atoms are found
    if all(iFound)
        m = iso.m(loc);
    end 
end

% aVogadro number (1/mol)
nA = 6.02214129e23;

% find formula
diffLabel = unique(aLabel);

numAtom = cell(1,2*numel(diffLabel));

for ii = 1:numel(diffLabel)
    temp = sum(strcmp(aLabel,diffLabel{ii}));
    numAtom{2*ii-1} = diffLabel{ii};
    numAtom{2*ii} = temp;
end

% number of formula units in the unit cell
nForm = double(obj.unit.nformula);
if nForm == 0
    % formula units is not given directly, try to find greates commond divider
    nForm = gcdv([numAtom{2:2:end}]);
end

% number of formula in cell
formula.N = nForm;
% formula mass in gramms
formula.m = sum(m)/nForm;
% cell volume in Angstrom^3
formula.V = det(obj.basisvector);
% density in g/cm^3
formula.rho = formula.m/formula.V/nA*1e24*nForm;

% divide the number of atoms in formula
numAtom(2:2:end) = num2cell([numAtom{2:2:end}]/nForm);

formula.chemform = sprintf('%s%d',numAtom{:});
formula.chemlabel = numAtom(1:2:end);
formula.chemnum = [numAtom{2:2:end}];

if nargout > 0
    varargout{1} = formula;
end

if nargout == 0
    fprintf0(fid,'     <strong>Chemical formula:</strong>  %s\n',formula.chemform);
    fprintf0(fid,'     <strong>Formula mass:</strong>      %8.3f g/mol\n',formula.m);
    fprintf0(fid,'     <strong>Formula in cell:</strong>   %8d units\n',formula.N);
    fprintf0(fid,'     <strong>Cell volume:</strong>       %8.3f Angstrom^3\n',formula.V);
    fprintf0(fid,'     <strong>Density:</strong>           %8.3f g/cm^3\n',formula.rho);
end

end