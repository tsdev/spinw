function varargout = powder_optimiser(obj,data,varargin)

%% Function workflow

% Remove the elastic line if it is there.

% We set all outside range to NaN.

% We chooose N q (NQ) slices to optimise on.

% Each Q slice will be of P points (QP).

% We start the pso optimiser with the current values.

% In each loop we will use the fibonaci option and have NFIB points.

% The computed sprectum will be convolved and compared to the data. R^2? 

% At the update stage `matparser` will be used to update the exchanges.

% When a solution is approximated a least squares method will be employed.

% The optimised object will be returned with an R value.


%% Input parameters.

% Initial spinW object

% Data to be optimised (single dat.x, dat.y, dat.z, dat.e for now)

% OPTIONAL

% Exchange labels to be fitted.

% Starting points.

% PSO options.

end

