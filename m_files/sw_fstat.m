function parOut = sw_fstat(state, parIn, T, E, M, ~)
% calculates termodynamical averages during an annealing simulation
%
% parOut = SW_FSTAT(state, parIn, T, E, M, nExt)
%
% It calculates statistical properties of different physical variables over
% several sample state. Called by sw.anneal.
%
% Input:
%
% state         Defines the task of the function.
%               1   Initialize the parOut structure.
%               2   Store the parameters of the physical state.
%               3   Calculate physical properties from the variable
%                   statistics.
% parIn         Same as parOut.
% T             Temperature of the system, vector [1,nT].
% E             Energy of the system, vector [1 nT].
% M             Magnetic moment of every atom, size: [spinDim,nMagExt*nT].
% nExt          Size of the magnetic cell, size: [3,1].
% kB            Boltmann constant, units of temperature.
%
% Output:
%
% parOut        Output parameter structure.
% parOut.nStat  The number of evaluated states.
% parOut.M      <M> summed over all magnetic moment, dimensions are
%               [spinDim,nMagExt*nT].
% parOut.M2     <M^2> summed over all magnetic moment, dimensions are
%               [spinDim,nMagExt*nT].
% parOut.E      <E> summed over all magnetic moment.
% parOut.E2     <E^2> summed over all magnetic moment.
%
%
% For the final execution, the following parameters are calculated:
%
% parOut        Array of struct, size [1 nT].
% parOut.avgM   Average components of the magnetisation over nStat runs,
%               size: (3,nMagExt).
% parOut.stdM   Standard deviation of the mgnetisation components over
%               nStat runs, size: (3,nMagExt).
% parOut.avgE   Average system energy per spin over nStat runs, scalar.
% parOut.stdE   Standard deviation of the system energy per spin over
%               nStat runs, scalar.
% parOut.T      Final temperature of the sample.
% parOut.Cp     Heat capacity of the sample: (<E^2>-<E>^2)/kB/T^2.
% parOut.Chi    Magnetic susceptibility of the sample: (<M^2>-<M>^2)/kB/T.
%
% See also SW, SW.ANNEAL.
%

if nargin == 0
    help sw_fstat;
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