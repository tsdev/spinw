function [formFactVal, coeff, S] = sw_mff(atomName, Q, nCoeff)
% returns the magnetic form factor values and coefficients
% 
% ### Syntax
% 
% `[~, coeff, s] = sw_mff(atomname)`
% 
% `[formfactval, coeff, s] = sw_mff(atomname,Q)`
%
% ### Description
% 
% `[~, coeff, s] = sw_mff(atomname)` returns the magnetic form
% factor coefficients for the magnetic atom identified by a string, e.g.
% `'MCR3'`. The function reads the `magion.dat` file for the stored form
% factor coefficients.
%
% `[formfactval, coeff, s] = sw_mff(atomname,Q)` also calculates the form
% factor values at the given $Q$ points (in \\Angstrom$^{-1}$ units.
%
% The source of the form factor data are:
% 1. A.-J. Dianoux and G. Lander, Neutron Data Booklet (2003).
% 2. K. Kobayashi, T. Nagao, and M. Ito, Acta Crystallogr. A. 67, 473 (2011).
%
% ### Input Arguments
% 
% `atomName`
% : String, contains the name of the magnetic ion in FullProf
%   notation (e.g. for Cr$^{3+} use `'MCR3'` or `'Cr3'`). It can be also a
%   vector of the 7 form factor coefficients. If the string contains
%   whitespace, the first word will be used as input. Can be also a cell of
%   strings to calculate coefficients for multiple ions.
% 
% `Q`
% : Momentum transfer in \\Angstrom$^{-1}$ units in a matrix with dimensions of
%   $[1\times n_Q]$ or $[3\times n_Q]$.
% 
% ### Output Arguments
% 
% `formFactVal`
% : Value of the form factor, evaluated at the given $Q$ points.
%
% `coeff`
% : Form factor coefficients according to the following formula:
%   
%   $\langle j_0(Q_s)\rangle = A\exp(-a\cdot Q_s^2) + B\exp(-b\cdot Q_s^2) + C\exp(-c\cdot Q_s^2) + D\exp(-d\cdot Q_s^2) + E$
%
%   where $Q_s = \frac{Q}{4\pi}$ and $A$, $a$, $B$, ... are the coefficients.
%   The $D$ and $d$ coefficients can be zero.
% 
% `S`
% : Value of the spin quantum number (read from the spin column in `magion.dat`).
%

% by default return 9 numbers as coefficients
nCoeff0 = 9;

if nargin == 0
    help sw_mff
    return
end

if nargin == 1
    Q = [];
end

if nargin < 3
    nCoeff = nCoeff0;
end



if ischar(atomName)
    atomName = {atomName};
end

if iscell(atomName)
    % open the form factor definition file
    ffPath = [sw_rootdir 'dat_files' filesep 'magion.dat'];
    % read form factor data
    formFact = sw_readtable(ffPath);
    % constant 1 form factor for atoms couldn't find
    formFact(end+1).a  = zeros(1,4);
    formFact(end).b    = zeros(1,4);
    formFact(end).c    = 1;
    formFact(end).spin = 0;
    
    idx = zeros(1,numel(atomName))+numel(formFact);
    
    for ii = 1:numel(atomName)
        
        % if there is whitespace, use the second word
        atomName0 = strword(atomName{ii},2,true);
        atomName0 = atomName0{1};
        
        % remove +/- symbols
        atomName0 = atomName0(atomName0>45);
        
        % search for the name of the atom
        idx0 = [find(strcmpi({formFact(:).label},atomName0))...
            find(strcmpi(cellfun(@(C)C(2:end),{formFact(:).label},'UniformOutput',0),atomName0))];
        
        if ~isempty(idx0)
            idx(ii) = idx0(1);
        end
    end
    
    coeff = [reshape([formFact(idx).a],4,[]);reshape([formFact(idx).b],4,[]);[formFact(idx).c]]';
    coeff = coeff(:,[1 5 2 6 3 7 4 8 9]);
    S     = [formFact(idx).spin];
    
    % pad the number of coefficients
    if nCoeff > nCoeff0
        coeff = [coeff(:,1:8) zeros(size(coeff,1),nCoeff-9) coeff(:,9)];
    end
    
    if any(idx == numel(formFact))
        fIdx = find(idx == numel(formFact));
        warning('sw_mff:WrongInput','The magnetic form factor for %s is undefined, constant 1 will be used instead!',atomName{fIdx(1)})
    end
elseif size(atomName,2) >= 9
    % just calculates the form factor values
    coeff = atomName;
else
    error('sw_mff:WrongInput','Wrong input!')
end

if ~isempty(Q)
    if size(Q,1)==3
        Q = sqrt(sum(Q.^2,1));
    elseif size(Q,1)~=1
        error('sw_mff:WrongInput','The dimensions of the given Q array is wrong!')
    end
    Qs = Q(:)'/(4*pi);
    
    formFactVal = zeros(size(coeff,1),numel(Qs));
    
    for ii = 2:2:size(coeff,2)
        formFactVal = formFactVal + bsxfun(@times,coeff(:,ii-1),exp(bsxfun(@times,-coeff(:,ii),Qs.^2)));
    end
    
    formFactVal = bsxfun(@plus,formFactVal,coeff(:,end));
    
else
    formFactVal = [];
end

end