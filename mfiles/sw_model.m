function obj = sw_model(model, param, fid)
% creates different predefined spin models
%
% obj = SW_MODEL(model, param, {fid})
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
%               'chain'     Chain with further neighbor interactions.
%
% param     Input parameters of the model, depending on which is selected.
% fid       Where to print the text output. Default is 1 to print to the
%           Command Window. Optional.
%
% Output:
%
% obj       spinw class object with the selected model.
%
% See also SPINW.
%

if nargin == 0
    help sw_model
    return
end

if nargin < 3
    fid = 1;
end

fprintf0(fid,'Preparing ''%s'' model ...\n',model);

obj = spinw;
fileid(obj,fid)

switch model
    case 'triAF'
        obj.genlattice('lat_const',[3 3 9],'angled',[90 90 120])
        obj.addatom('r',[0 0 0],'S',1,'color','darkmagenta')
        obj.gencoupling('maxDistance',10)
        if nargin > 1
            for ii = 1:numel(param)
                obj.addmatrix('value',param(ii),'label',['J' num2str(ii)],'color',sw_colorname(randi(140),1))
                obj.addcoupling('mat',ii,'bond',ii)
            end
            
            %obj.lattice.lat_const(3) = 5;
            obj.genmagstr('mode','helical','S',[1 0 0]','k',[0 0 0])
            obj.optmagstr('func',@gm_planar,'xmin',[0 0 0 0 0 0],'xmax',[0 1/2 1/2 0 0 0],'nRun',10)
        end
    case 'squareAF'
        obj.genlattice('lat_const',[3 3 9],'angled',[90 90 90])
        obj.addatom('r',[0 0 0],'S',1,'color','darkmagenta')
        obj.gencoupling('maxDistance',10)
        if nargin>1
            for ii = 1:numel(param)
                obj.addmatrix('value',param(ii),'label',['J' num2str(ii)],'color',sw_colorname(randi(140),1))
                obj.addcoupling('mat',ii,'bond',ii)
            end
            
            %obj.lattice.lat_const(3) = 5;
            obj.genmagstr('mode','helical','S',[1 0 0]','k',[0 0 0])
            obj.optmagstr('func',@gm_planar,'xmin',[0 0 0 0 0 0],'xmax',[0 1/2 1/2 0 0 0],'nRun',10)
            tol = 2e-4;
            helical = sum(abs(mod(abs(2*obj.magstr.k)+tol,1)-tol).^2) > tol;
            if ~helical
                obj.genmagstr('mode','helical','nExt',[2 2 1]);
            end
        end
    case 'chain'
        obj.genlattice('lat_const',[3 9 9],'angled',[90 90 90])
        obj.addatom('r',[0 0 0],'S',1,'color','darkmagenta')
        obj.gencoupling('maxDistance',10)
        if nargin>1
            for ii = 1:numel(param)
                obj.addmatrix('value',param(ii),'label',['J' num2str(ii)],'color',sw_colorname(randi(140),1))
                obj.addcoupling('mat',ii,'bond',ii)
            end
            
            %obj.lattice.lat_const(2:3) = 5;
            obj.genmagstr('mode','helical','S',[1 0 0]','k',[0 0 0])
            obj.optmagstr('func',@gm_planar,'xmin',[0 0 0 0 0 0],'xmax',[0 1/2 0 0 0 0],'nRun',10)
            tol = 2e-4;
            helical = sum(abs(mod(abs(2*obj.magstr.k)+tol,1)-tol).^2) > tol;
            if ~helical && any(obj.magstr.k>tol)
                nExt = [1 1 1];
                nExt(obj.magstr.k(1:2)>tol) = 2;
                obj.genmagstr('mode','helical','nExt',nExt);
            end
        end
    otherwise
        error('sw_model:WrongINput','Model does not exists!')
end

fprintf0(fid,'... ready!\n');

end