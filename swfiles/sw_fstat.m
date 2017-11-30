function parOut = sw_fstat(state, parIn, T, E, M, ~)
% calculates thermodynamical averages
% 
% ### Syntax
% 
% `parOut = sw_fstat(state, parIn, T, E, M, nExt)`
% 
% ### Description
% 
% `parOut = sw_fstat(state, parIn, T, E, M, nExt)` calculates statistical
% properties of different physical variables over several sampled state.
% The function is called by [spinw.anneal].
% 
% ### Input Arguments
% 
% `state`
% : Defines the task of the function.
%   * `1`   Initialize the parOut structure.
%   * `2`   Store the parameters of the physical state.
%   * `3`   Calculate physical properties from the variable
%           statistics.
% 
% `parIn`
% : Same as `parOut`.
% 
% `T`
% : Temperature of the system, row vector with $n_T$ number of elements.
% 
% `E`
% : Energy of the system, row vector with $n_T$ number of elements.
% 
% `M`
% : Magnetic moment of every atom in a matrix with dimensions of $[d_{spin}\times n_{magExt}\cdot n_T]$.
% 
% `nExt`
% : Size of the magnetic supercell, column vector of 3 integers.
% 
% `kB`
% : Boltzmann constant, units of temperature.
% 
% ### Output Arguments
% 
% `parOut`
% : Output parameter structure with the following fields:
%   * `nStat`   The number of evaluated states.
%   * `M`       $\langle M\rangle$ averaged over all magnetic moment stored
%               in a matrix with dimensions of $[d){spin}\times
%               n_{magExt}\cdot n_T]$.
%   * `M2`      $\langle M^2\rangle$ averaged over all magnetic moment
%               stored in a matrix with dimensions of $[d){spin}\times
%               n_{magExt}\cdot n_T]$.
%   * `E`       $\langle E\rangle$  summed over all magnetic moment.
%   * `E2`      $\langle E^2\rangle$  summed over all magnetic moment.
%
% For the final execution, the following parameters are calculated:
% `parOut`
% : Array of struct with $n_T$ number of elements:
%   * `avgM`    Average components of the magnetisation over $n_{stat}$ runs,
%               matrix with dimensions of $[3\times n_{magExt}]$.
%   * `stdM`    Standard deviation of the mgnetisation components over
%               $n_{stat}$ runs, matrix with dimensions of $[3\times n_{magExt}]$.
%   * `avgE`    Average system energy per spin over $n_{stat}$ runs, scalar.
%   * `stdE`    Standard deviation of the system energy per spin over
%               $n_{stat}$ runs, scalar.
%   * `T`       Final temperature of the sample.
%   * `Cp`      Heat capacity of the sample: $(\langle E^2\rangle-\langle E\rangle^2)/k_B/T^2$.
%   * `Chi`     Magnetic susceptibility of the sample: $(\langle M^2\rangle-\langle M\rangle^2)/k_B/T$.
% 
% ### See Also
% 
% [spinw.anneal]
%

if nargin == 0
    swhelp sw_fstat;
    return
end

switch state
    case 2
        % Save parameter statistics.
        parOut.nStat = parIn.nStat + 1;
        parOut.M     = parIn.M  + M;
        parOut.M2    = parIn.M2 + M.^2;
        parOut.E     = parIn.E  + E;
        parOut.E2    = parIn.E2 + E.^2;
    case 1
        % Initialise parameters.
        parOut.nStat = 0;
        parOut.M     = 0*M;
        parOut.M2    = 0*M;
        parOut.E     = 0*T;
        parOut.E2    = 0*T;
    case 3
        % Calculates the mean and standard deviation of magnetic moment and energy.
        spinDim    = size(M,1);
        nMagExt    = size(M,2)-1;
        nStat      = parIn.nStat;
        kB         = parIn.kB;
        
        parIn.M2   = parIn.M2(:,1:end-1);
        parIn.M    = parIn.M(:,1:end-1);
        % CHECK
        parIn.M2(spinDim+1:3,:) = 0;
        parIn.M(spinDim+1:3,:)  = 0;
        
        % Save all statistical data to the temporary parTemp struct.
        parOut.avgM = parIn.M/nStat;
        parOut.stdM = sqrt((parIn.M2 - parIn.M.^2/nStat)/nStat);
        parOut.avgE = parIn.E/nStat/nMagExt;
        parOut.stdE = sqrt((parIn.E2/nMagExt - parIn.E.^2/nMagExt/nStat)/nMagExt/nStat);
        parOut.T    = T;
        
        parOut.M    = M(:,1:end-1);
        parOut.M(spinDim+1:3,:) = 0;
        
        % Heat capacity per magnetic moment.
        parOut.Cp   = parOut.stdE.^2./kB./(T.^2);
        parOut.Chi  = sum(parOut.stdM.^2/nMagExt,2)/kB/T;
        
end
end