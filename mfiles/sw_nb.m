function bc = sw_nb(atomName)
% returns the bound coherent neutron scattering length (fm)
%
% bc = SW_NB(atomName)
%
% Input:
%
% atomName      String, contains the name of the atom or isotope (e.g. 13C for the isotope of carbon).
%
% Output:
%
% bc        Value of the bound coherent neutron scattering length  in fm.
%

if nargin == 0
    help sw_nb
    return
end

if ischar(atomName)
    atomName = {atomName};
end

if iscell(atomName)
    % open the form factor definition file
    ffPath = [sw_rootdir 'dat_files' filesep 'isotope.dat'];
    % read neutron scattering length data
    sigma = sw_readtable(ffPath);
    % constant 1 scattering length for atoms couldn't find
    sigma(end+1).bc = 1;
    sigma(end).name = 'Un';
    
    idx = zeros(1,numel(atomName))+numel(sigma);
    
    for ii = 1:numel(atomName)
        
        % if there is whitespace, use the second word
        atomName0 = strword(atomName{ii},2,true);
        atomName0 = atomName0{1};
        
        % remove +/- symbols
        atomName0 = atomName0(atomName0>45);
        
        % extract atomic mass value if given
        A0 = sscanf(atomName0,'%d');
        if isempty(A0)
            A0 = -1;
        end
        
        % remove all numbers
        atomName0 = atomName0(atomName0>'9');
        
        % remove M for magnetic atoms
        if atomName0(1) == 'M' && numel(atomName0)>1 && upper(atomName0(2))==atomName0(2)
            atomName0 = atomName0(2:end);
        end
        % search for the name of the atom
        idx0 = find(strcmpi({sigma(:).name},atomName0));
        
        % search the right mass
        idx0 = idx0([sigma(idx0).A]==A0);
        
        if ~isempty(idx0)
            idx(ii) = idx0;
        end
    end
    
    bc = [sigma(idx).bc];
    
    if any(idx == numel(sigma))
        fIdx = find(idx == numel(sigma));
        warning('sw_nb:WrongInput','The scattering length for %s is undefined, constant 1 will be used instead!',atomName{fIdx(1)})
    end
else
    error('sw_mff:WrongInput','Wrong input!')
end

end