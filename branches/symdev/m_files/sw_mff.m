function [formFactVal, coeff, S] = sw_mff(atomName, Q)
% [formFactVal, coeff, S] = SW_MFF(atomName, {Q}) returns the magnetic
%  form factor values and the coefficients.
%
% Input:
%
% atomName      String, contains the name of the magnetic ion in FullProf
%               notation (e.g. Cr^3+ --> 'MCR3' or 'Cr3'). It can be also a
%               vector of the 7 coefficients, see below.
% Q             Momentum transfer in Angstrom^-1 units, optional.
%
% Output:
%
% formFactVal   Value of form factor, evaluated at the Q points if Q is
%               defined.
% coeff         Form factor coefficients according to the following
%               formula:
%               <j0(Qs)> = A*exp(-a*Qs^2) + B*exp(-b*Qs^2) + C*exp(-c*Qs^2) + D,
%               where Qs = Q/(4*pi) and A, a, B, ... are the coefficients.
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
    while ~feof(fid)
        fLine = fgets(fid);
        [formFact.name{idx,1}, ~, ~, nextIdx] = sscanf(fLine,'%s',1);
        [formFact.coeff(idx,1:7), ~, ~, nextIdx2] = sscanf(fLine(nextIdx:end),'%f,',7);
        formFact.S(idx) = sscanf(fLine(nextIdx+nextIdx2:end),'%f,',1);
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
        S     = formFact.S(idx);
    else
        coeff = [zeros(1,6) 1];
    end
    
else
    %flag = true;
    coeff = atomName;
end



if nargin > 1
    Qs = Q/(4*pi);
    formFactVal = coeff(1)*exp(-coeff(2)*Qs.^2)+coeff(3)*exp(-coeff(4)*Qs.^2)+...
        coeff(5)*exp(-coeff(6)*Qs.^2)+coeff(7);
else
    formFactVal = [];
end

end