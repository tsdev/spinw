function obj = sw_model(model, param)
% obj = SW_MODEL(model, param) creates different predefined spin models.
%
% Input:
%
% model     String, name of the model, one of the following:
%               'triAF'     Triangular lattice Heisenberg antiferromagnet
%                           in the ab plane (a=b=3 Angstrom), with gamma =
%                           120 deg angle and optimised magnetic structure.
%                           Arbitrary number of Heisenberg interaction can
%                           be defined, param(1) gives the value of 1st
%                           neighbor interaction, param(2) the second etc.
%               'squareAF'  Square lattice antiferromagnet.
%               'chainAF'   Antiferromagnetic chain.
%
% param     Input parameters of the model, depending on which is selected.
%
% Output:
%
% obj       sw class object with the selected model.
%
% See also SW.
%

fprintf('Preparing ''%s'' model ...',model)

switch model
    case 'triAF'
        obj = sw;
        obj.genlattice('lat_const',[3 3 500],'angled',[90 90 120])
        obj.addatom('r',[0 0 0],'S',1,'color','darkmagenta')
        obj.gencoupling('maxDistance',10)
        
        for ii = 1:numel(param)
            obj.addmatrix('value',param(ii),'label',['J' num2str(ii)],'color',sw_colorname(randi(140),1))
            obj.addcoupling(ii,ii)
        end
        
        obj.lattice.lat_const(3) = 5;
        obj.genmagstr('mode','direct','S',[1 0 0])
        obj.optmagstr('func',@gm_planar,'xmin',[0 0 0 0 0 0],'xmax',[0 1/2 1/2 0 0 0],'nRun',10,'fid',0)
    case 'squareAF'
        obj = sw;
        obj.genlattice('lat_const',[3 3 500],'angled',[90 90 90])
        obj.addatom('r',[0 0 0],'S',1,'color','darkmagenta')
        obj.gencoupling('maxDistance',10)
        
        for ii = 1:numel(param)
            obj.addmatrix('value',param(ii),'label',['J' num2str(ii)],'color',sw_colorname(randi(140),1))
            obj.addcoupling(ii,ii)
        end
        
        obj.lattice.lat_const(3) = 5;
        obj.genmagstr('mode','direct','S',[1 0 0])
        obj.optmagstr('func',@gm_planar,'xmin',[0 0 0 0 0 0],'xmax',[0 1/2 1/2 0 0 0],'nRun',10,'fid',0)
        
    otherwise
        error('sw_model:WrongINput','Model does not exists!')
end

fprintf(' ready!\n')

end