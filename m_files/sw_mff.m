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
%               Optionally the (D,d) coefficients are zero.
%
% S             Value of the spin quantum number (last column in ion.dat).
%

S = 0;

if nargin == 0
    help sw_mff;
    return;
end

formFact = struct;

if ischar(atomName)
    % if there is whitespace, use the first word
    atomName = strword(atomName,1);
    atomName = atomName{1};
    %     atomName = strsplit(atomName);
    %     wordIdx = find(cellfun(@numel,atomName));
    %     atomName = atomName{wordIdx(1)};
    
    % remove +/- symbols
    atomName = atomName(atomName>45);
    % Open the form factor definition file.
    ffPath = [sw_rootdir 'dat_files' filesep 'ion.dat'];
    fid = fopen(ffPath);
    if fid == -1
        error('sw:sw_mff:FileNotFound',['Form factor definition file not found: '...
            regexprep(ffPath,'\' , '\\\') '!']);
    end
    
    % Read all line from file
    idx = 1;
    formFact.coeff = zeros(0,15);
    
    while ~feof(fid)
        fLine = fgets(fid);
        [formFact.name{idx,1}, ~, ~, nextIdx] = sscanf(fLine,'%s',1);
        [coeffTemp, ~, ~, nextIdx2] = sscanf(fLine(nextIdx:end),'%f ',inf);
        formFact.coeff(idx,1:numel(coeffTemp)-1) = coeffTemp(1:end-1);
        formFact.S(idx) = coeffTemp(end);
        idx = idx + 1;
    end
    fclose(fid);
    
    flag = false;
    
    % Looks for the name of the atom.
    for ii = 1:length(formFact.name)
        selected = formFact.name{ii};
        if ~isempty(strfind(upper(atomName),upper(selected))) || ~isempty(strfind(upper(atomName),upper(selected(2:end))))
            idx = ii;
            flag = true;
        end
    end
    
    if flag
        coeff = formFact.coeff(idx,:);
        % remove trailing zeros
        cutIdx = find([diff(coeff==0) 0],1,'last');
        coeff = coeff(1:cutIdx);
        
        S     = formFact.S(idx);
    else
        coeff = [zeros(1,6) 1];
    end
    
else
    %flag = true;
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

end