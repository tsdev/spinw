function [formFactVal, coeff, label] = sw_cff(atomName, Q)
% returns the atomic charge form factor values for X-ray
%
% [formFactVal, coeff, label] = SW_CFF(atomName, {Q})
%
% The provided form factor values at Q=0 are normalized to Z.
%
% Input:
%
% atomName      String, contains the name of the ion in using the symbol of
%               the element following the charge, e.g. Cr3+. It can be also
%               the coefficients to calculate f. If the string contains
%               whitespace, the first word will be used as input.
% Q             Momentum transfer in Angstrom^-1 units with dimensions of
%               [1 nQ] or [3 nQ], optional.
%
% Output:
%
% formFactVal   Value of form factor, evaluated at the Q points if Q is
%               defined.
% coeff         Form factor coefficients according to the following
%               formula:
%                   f0(Qs) = c + SUM a_i*EXP(-b_i*(Qs^2))
%                                i=1,5
%               where Qs = Q/(4*pi) and [a_1 b_1 a_2 b_2 ... c] are the
%               coefficients.
%

if nargin == 0
    help sw_cff
    return
end

if ischar(atomName)
    % if there is whitespace, use the first word
    atomName = strword(atomName,1);
    atomName = atomName{1};
    
    % open the form factor definition file
    ffPath = [sw_rootdir 'dat_files' filesep 'xrayion.dat'];
    % read form factor data
    formFact = sw_readtable(ffPath);
    
    flag = false;
    
    % search for the name of the atom
    for ii = 1:numel(formFact)
        selected = formFact(ii).label;
        if strcmpi(atomName,selected)
            idx = ii;
            flag = true;
        end
    end
    
    if flag
        coeff = [formFact(idx).a; formFact(idx).b];
        coeff = [coeff(:)' formFact(idx).c];
        if isempty(coeff)
            coeff = 0;
        end
        
    else
        coeff = [zeros(1,11) 1];
        if nargout < 3
            warning('sw_cff:WrongInput',['The form factor for %s is '...
                'undefined, constant 1 will be used instead!'],atomName)
        end
    end
    
else
    coeff = atomName;
end

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

if nargout > 2
    if flag
        label = formFact(idx).label;
    else
        label = [];
    end
end

end