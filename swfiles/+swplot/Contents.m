% package for 3D plotting
%
% The package contains functions that can create and control plotting of 3D
% objects. It contains low level functions to plot specific shapes on any
% figure (cylinder, line, arrow, ellipsoid, polyhedron, sphere, text) that
% are vectorized for fast plotting of multiple objects using a single
% [matlab:patch] command. Moreover it provides the [swplot.plot] command to
% plot multiple objects with specific color, transparency, label and legend
% on the swplot figure including tha ability for smooth rotation and
% translations via the mouse (that are much better than the default 3D
% rotation with mouse in Matlab). Moreover it has high level commands to
% plot crystals from SpinW objects (`swplot.plot...` commands).
%
% ### Files
%
% #### Basic 3D objects
%
% Plots primitive geometrical shapes, able to plot multiple objects when
% input is vectorized. Can give considerable speedup because a single patch
% command is used to plot multiple shapes.
%
%   swplot.arrow     
%   swplot.circle    
%   swplot.cylinder  
%   swplot.ellipsoid 
%   swplot.line      
%   swplot.polyhedron
%   swplot.text      
%
% #### Creating and modifying swplot figure
%
% These functions provide complete control over the transformations,
% objects on an swplot figure. These are the alternatives to the Matlab
% built in commands that control 3D plot axes. The advantage here is the
% smooth mouse rotation/zoom functionality, nice legend and ability to
% assign a tooltip text to any object that is triggered by a mouse click on
% the object.
%
%   swplot.activefigure
%   swplot.add         
%   swplot.base        
%   swplot.clear       
%   swplot.close       
%   swplot.delete      
%   swplot.export      
%   swplot.figure      
%   swplot.findobj     
%   swplot.getdata     
%   swplot.legend      
%   swplot.mouse       
%   swplot.plot        
%   swplot.tooltip     
%   swplot.transform   
%   swplot.translate   
%   swplot.view        
%   swplot.zoom        
%
% #### Plotting SpinW object on swplot figure
%
% These high level plot commands are used to plot different types of data
% that is stored in the [spinw] object. All inputs of these commands can be
% controlled from the [spinw.plot] command.
%
%   swplot.plotatom
%   swplot.plotbase
%   swplot.plotbond
%   swplot.plotcell
%   swplot.plotchem
%   swplot.plotion 
%   swplot.plotmag 
%
% #### Related functions
%
% Helper function that can be also used as standalone.
%
%   swplot.color          
%   swplot.icomesh        
%   swplot.logo           
%   swplot.subfigure
%   swplot.subplot
%