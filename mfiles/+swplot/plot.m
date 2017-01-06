function plot(varargin)
% plots objects to swplot figure
%
% SWPLOT.PLOT()
%
% Options:
%
% type      Type of object to plot in a string. Possible options are:
%               'arrow'         position specifies start and end points
%               'ellipsoid'     position specifies center
%               'cylinder'      position specifies start and end points
%               'circle'        position specifies center and normal vector
%               'line'          position specifies start and end points
% position  Position of the object/objects in a matrix with dimensions of
%           [3 2 nObject]/[3 1 nObject] depending on the type of object. 
% text      Text to appear in legend in a string for the same text of all 
%           objects or strings in a cell for multiple objects with
%           dimension [1 nObject]. Default is empty string.
% legend    If true, the objects appear in the legend single as a single
%           entry (if all objects has the same text and color) or multiple
%           entry (otherwise). Default is false.
% color     Color of objects, either a single color or as many colors as
%           many objects are given in a matrix with dimensions of [3 1]/[3
%           nObject]. Default is red.
% units     
% 


inpForm.fname  = {'type' 'position' 'text' 'legend' 'color' 'units' 'name'};
inpForm.defval = {false     false    true       0        1e-4  true    };
inpForm.size   = {[1 1]     [1 1]    [1 1]      [1 1]    [1 1] [1 1]   };

% inpForm.fname  = [inpForm.fname  {'omega_tol' 'saveSabp' 'saveV' 'saveH'}];
% inpForm.defval = [inpForm.defval {1e-5        true       false   false  }];
% inpForm.size   = [inpForm.size   {[1 1]       [1 1]      [1 1]   [1 1]  }];

param = sw_readparam(inpForm, varargin{:});




end