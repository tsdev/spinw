% general functions
%
% This folder contains all the spectral functions and general functions
% that are related to SpinW.
%
% ### Files
%
% #### Transforming and plotting calculated spin wave spectrum
%
% These functions operate on the calculated spectra, which is the output of
% [spinw.spinwave] or [spinw.powspec] commands. They enable to post process
% the calculated spin-spin correlation function, including instrumental
% resolution, cross section calculation, binning etc.
%
%   sw_econtract 
%   sw_egrid     
%   sw_filelist  
%   sw_instrument
%   sw_magdomain 
%   sw_neutron   
%   sw_omegasum  
%   sw_plotspec  
%   sw_xray      
%
% #### Generate list of vectors in reciprocal space
%
% These two functions can generate a set of 3D points in reciprocal space
% defining either a path made out of straigh lines or a volume.
%
%   sw_qgrid
%   sw_qscan
%
% #### Resolution claculation and convolution
%
% These functions can import Energy resolution function and convolute it
% with arbitrary multidimensional dataset
%
%   sw_res    
%   sw_resconv
%   sw_tofres 
%
% #### SpinW model related functions
%
%   sw_extendlattice
%   sw_fstat        
%   sw_model        
%   sw_bonddim      
%
% #### Constraint functions
%
% Contraint functions for [spinw.optmagstr].
%
%   gm_planar      
%   gm_planard     
%   gm_spherical3d 
%   gm_spherical3dd
%
% #### Geometrical calculations
%
% Basic geometrical calculators, functions to generatate rotation
% operators, generate Cartesian coordinate system from a set of vectors,
% calculate normal vector to a set of vector, etc.
%
%   sw_basismat
%   sw_cartesian
%   sw_fsub     
%   sw_mattype  
%   sw_nvect    
%   sw_quadell  
%   sw_mirror
%   sw_rot      
%   sw_rotmat   
%   sw_rotmatd  
%
% #### Text and graphical input/output for different high level commands
%
%   sw_multicolor  
%   sw_parstr      
%   sw_timeit      
%
% #### Acessing the SpinW database
%
% Functions to read the different data files that store information on
% atomic properties, such as magnetic form factor, charge, etc.
% 
%   sw_atomdata
%   sw_cff     
%   sw_mff     
%   sw_nb      
%
% #### Useful physics functions
%
% The two functions can calculate the Bose factor and convert
% energy/momentum units, both usefull for neutron and x-ray scattering.
%
%   sw_bose     
%   sw_converter
%
% #### Import functions
%
% Functions to import tables in text format.
%
%   sw_import   
%   sw_readspec 
%   sw_readtable
%
% #### Miscellaneous
%
%   swdoc
%   sw_freemem   
%   sw_readparam 
%   sw_rootdir   
%   sw_uniquetol 
%   sw_update    
%   sw_version   
%   sw_mex       
