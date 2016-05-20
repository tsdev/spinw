function [formFactVal, coeff, S] = sw_mff(atomName, Q)
% returns the magnetic form factor values and the coefficients
%
% [formFactVal, coeff, S] = SW_MFF(atomName, {Q})
%
% Input:
%
% atomName      String, contains the name of the magnetic ion in FullProf
%               notation (e.g. Cr^3+ --> 'MCR3' or 'Cr3'). It can be also a
%               vector of the 7 coefficients, see below. If the string
%               contains whitespace, the first word will be used as input.
% Q             Momentum transfer in Angstrom^-1 units with dimensions of
%               [1 nQ] or [3 nQ], optional.
%
% Output:
%
% formFactVal   Value of form factor, evaluated at the Q points if Q is
%               defined.
% coeff         Form factor coefficients according to the following
%               formula:
%               <j0(Qs)> = A*exp(-a*Qs^2) + B*exp(-b*Qs^2) + C*exp(-c*Qs^2) + D*exp(-d*Qs^2) + E,
%               where Qs = Q/(4*pi) and A, a, B, ... are the coefficients.
%               The (D,d) coefficients can be zero.
%
% S             Value of the spin quantum number (last column in magion.dat).
%
% The source for the form factor data are:
% [1] A.-J. Dianoux and G. Lander, Neutron Data Booklet (2003).
% [2] K. Kobayashi, T. Nagao, and M. Ito, Acta Crystallogr. A. 67, 473 (2011).
%

S = 0;

if nargin == 0
    help sw_mff
    return
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
    formFact(end+1).a = zeros(1,4);
    formFact(end).b   = zeros(1,4);
    formFact(end).c   = 1;
    
    idx = zeros(1,numel(atomName))+numel(formFact);
    
    for ii = 1:numel(atomName)
        
        % if there is whitespace, use the first word
        atomName0 = strword(atomName{ii},1);
        atomName0 = atomName0{1};
        
        % remove +/- symbols
        atomName0 = atomName0(atomName0>45);
        
        % search for the name of the atom
        idx0 = [find(strcmpi({formFact(:).label},atomName0))...
            find(strcmpi(cellfun(@(C)C(2:end),{formFact(:).label},'UniformOutput',0),atomName0))];
        
        if ~isempty(idx0)
            idx(ii) = idx0;
        end
    end
    
    coeff = [reshape([formFact(idx).a],4,[]);reshape([formFact(idx).b],4,[]);[formFact(idx).c]]';
    
    if nargout < 3 && any(idx == numel(formFact))
        fIdx = find(idx == numel(formFact));
        warning('sw_mff:WrongInput','The form factor for %s is undefined, constant 1 will be used instead!',atomName{fIdx(1)})
    end
    
else
    coeff = atomName;
end

% TODO for multiple atoms
if nargin > 1
    if all(size(Q)>1)
        % if Q points are given as a list of Q vectors in Angstrom^-1
        % units, calculate the absolute value of Q
        Q = sqrt(sum(Q.^2,1));
    end
    Qs = Q/(4*pi);
    
    formFactVal = Qs*0;
    for ii = 2:2:numel(coeff)
        formFactVal = formFactVal + coeff(ii-1)*exp(-coeff(ii)*Qs.^2);
    end
    formFactVal = formFactVal + coeff(end);
else
    formFactVal = [];
end

end