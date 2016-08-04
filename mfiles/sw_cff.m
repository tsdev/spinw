function [formFactVal, coeff] = sw_cff(atomName, Q)
% returns the atomic charge form factor values for X-ray scattering
%
% [formFactVal, coeff] = SW_CFF(atomName, {Q})
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
% See also sw_mff.
%

if nargin == 0
    help sw_cff
    return
end

if ischar(atomName)
    atomName = {atomName};
end

if iscell(atomName)
    % open the form factor definition file
    ffPath = [sw_rootdir 'dat_files' filesep 'xrayion.dat'];
    % read form factor data
    formFact = sw_readtable(ffPath);
    
    % constant 1 form factor for atoms couldn't find
    formFact(end+1).a  = zeros(1,5);
    formFact(end).b    = zeros(1,5);
    formFact(end).c    = 1;
    formFact(end).spin = 0;
    
    idx = zeros(1,numel(atomName))+numel(formFact);
    
    for ii = 1:numel(atomName)
        
        % if there is whitespace, use the second word
        atomName0 = strword(atomName{ii},2,true);
        atomName0 = atomName0{1};
        
        % remove leading uppercase letters up to the last one
        % TODO
        cutIdx = find(atomName0>='A' & atomName0<='Z',1,'last');
        if ~isempty(cutIdx) && cutIdx>1
            atomName0 = atomName0(cutIdx:end);
        end
        % add + symbol if not given for oxidation states
        if ~any(ismember(atomName0,'+-')) && any(ismember(atomName0,'1':'9'))
            atomName0 = [atomName0 '+']; %#ok<AGROW>
        end
        
        % search for the name of the atom
        idx0 = find(strcmpi({formFact(:).label},atomName0));
        
        atomName{ii} = atomName0;
        
        if ~isempty(idx0)
            idx(ii) = idx0;
        end
    end
    
    coeff = [reshape([formFact(idx).a],5,[]);reshape([formFact(idx).b],5,[]);[formFact(idx).c]]';
    coeff = coeff(:,[1 6 2 7 3 8 4 9 5 10 11]);
    
    
    if any(idx == numel(formFact))
        fIdx = find(idx == numel(formFact));
        warning('sw_cff:WrongInput','The x-ray scattering form factor for %s is undefined, constant 1 will be used instead!',atomName{fIdx(1)})
    end
elseif size(atomName,2) == 11
    coeff = atomName;
else
    error('sw_cff:WrongInput','Wrong input!')
end


if nargin > 1
    if all(size(Q)>1)
        % if Q points are given as a list of Q vectors in Angstrom^-1
        % units, calculate the absolute value of Q
        Q = sqrt(sum(Q.^2,1));
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