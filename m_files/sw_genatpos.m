function [rSave, symName] = sw_genatpos(sym, r, varargin)
% [rSave, symName] = SW_GENATPOS(sym, r, {print}) generates all symmetry 
% equivalent atomic positions from a given symmetry number and coordinates 
% of the first atom (r). If print is defined, the result is printed onto 
% the command window.
%
% Input:
%
% sym           Line index in the symmetry.dat file or string of the
%               symmetry operators.
% r             Atomic position in lattice units, dimensions are [1 3] or 
%               [3 1].
% {print}       Optional input, if defined, the result will be plotted on
%               the command window.
%
% Output:
%
% rSave         Generated atomic positions, dimensions are [3 nAtom].
% symName       String, the name of  the space group.
%
% See also SW, SW_GENCOUPLING.
%

[transf, transl, symName] = sw_gensym(sym);

nOp   = size(transf,3);
rSave = zeros(3,1);

% First atom is the original atom.
rSave(:,1) = mod(r,1);

Counter = 1;
while Counter<=size(rSave,2)

    
    for ii = 1:nOp 
        rNew  = mod(squeeze(transf(:,:,ii))*rSave(:,Counter)+transl(:,ii),1);        
        isnew = 1;
        for jj=1:size(rSave,2)
            %if sum(abs(Rnew-rSave(:,jj)))<1e-5
            % TODO fix this line
            if sum(1-2*abs(abs(rNew-rSave(:,jj))-0.5)) < 0.05
                isnew = 0;
            end
        end
        if isnew
           rSave(:,end+1) = rNew; %#ok<AGROW>
        end                   
    end
    Counter = Counter+1;
end

if nargin>2
    fprintf('Symmetry: %s\n',symName);
    fprintf('Atomic coordinates generated for: (%5.3f %5.3f %5.3f)\n',r);
    for ii = 1:size(rSave,2)
        fprintf('r%i (%5.3f %5.3f %5.3f)\n',ii,rSave(:,ii));
    end
end