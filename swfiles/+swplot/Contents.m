% package for 3D plotting
%
% The package contains functions that can create and control plot of 3D
% objects. It contains low level function to plot specific shapes on any
% figure (cylinder, line, arrow, ellipsoid, text) that are vectorized.
% Moreover it provides the plot command to plot multiple objects with
% specified color, transparency, label and legend on the swplot figure that
% has rotation and translations built in that are better than the default
% 3D rotation with mouse in Matlab. Moreover it has high level commands to
% plot crystals from SpinW objects (plotatom, plotmag, plotion, plotbond,
% plotbase, plotcell, plotchem).
%
% ### Files
%
% #### Basic 3D objects
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
%   swplot.ishg        
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
%   swplot.color          
%   swplot.icomesh        
%   swplot.logo           
%   swplot.patchfacefcn   
%   swplot.raytriangle    
%   swplot.setrangegui    
%   swplot.subfigure      
%   swplot.tooltipcallback
%   swplot.tooltipstring  
%