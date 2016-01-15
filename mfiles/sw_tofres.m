function spec = sw_tofres(spec, varargin)
% includes Q resolution to the spectrum
%
% spectra = SW_TOFRES(spectra, 'Option1', Value1, ...)
%
% Simulates the finite bin size of the cuts of TOF data.
%
% Input:
%
% spectra   Input structure, contains calculated correlation functions
%           withouth the resolution effect.
%
% Options:
%
% method    String, determines the method to genera the Q points, options:
%               'random'    The bin volume will be randomly sampled.
%               'grid'      The bin volume will be split up to a regular
%                           grid.
% dQ        Vector with three numbers of scalar. The width of the Q bin
%           along the three reciprocal lattice directions. The spectrum
%           will be integrated in the Q+/-(dQ/2) range. DEfault value is
%           [0.1 0.1 0.1].
% nQ        Vector with three numbers or scalar. Gives the number of Q
%           points along the three reciprocal lattice directions to average
%           over or the number of random Q points for the random method.
%
%
% Output:
%
% spectra that contains the calculated intensity in the swConv field.
%
% See also SW_EGRID, SW_INSTRUMENT.
%

% help when executed without argument
if nargin==0
    help sw_ressim
    return
end

dQ0 = ones(1,3)*0.1;
nQ0 = ones(1,3)*5;

inpForm.fname  = {'method'  'dQ'  'nQ'};
inpForm.defval = {'grid'    dQ0    nQ0};
inpForm.size   = {[1 -1] [1 -2] [1 -3]};


param = sw_readparam(inpForm, varargin{:});
obj   = spec.obj;
fid   = obj.fileid;

hkl   = spec.hkl;
Evect = spec.Evect;
nE    = numel(Evect)-1;
dQ    = param.dQ;
nQ    = param.nQ;

switch param.method
    case 'random'
        
        if numel(param.dQ) == 1
            dQ = ones(1,3)*dQ;
        end
        
        conv0 = zeros(nE,size(hkl,2));
        
        obj.fileid(0);
        
        fprintf0(fid,'Calculating TOF Q-binning using random method...\n')
        
        sw_status(0,1,fid);
        
        for ii = 1:prod(nQ)
            hklC = bsxfun(@minus,bsxfun(@times,dQ',rand(size(hkl))),dQ'/2)+hkl;
            spec0 = obj.spinwave(hklC,'formfact',spec.formfact);
            spec0 = sw_egrid(spec0,'Evect',Evect,'component',spec.component);
            conv0 = conv0 + spec0.swConv;
            
            sw_status(ii/prod(nQ)*100,0,fid);
        end
        
    case 'grid'
        
        
        if numel(param.dQ) == 1
            dQ = ones(1,3)*dQ;
        end
        
        if numel(nQ) == 1
            nQ = ones(1,3)*nQ;
        end
        
        % vector along each rl direction
        Qv = cell(1,3);
        for ii = 1:3
            Qv{ii}= linspace(-dQ(ii),dQ(ii),nQ(ii));
        end
        
        dhkl = cell(1,3);
        [dhkl{1},dhkl{2},dhkl{3}] = ndgrid(Qv{:});
        
        
        conv0 = zeros(nE,size(hkl,2));
        
        obj.fileid(0);
        
        fprintf0(fid,'Calculating TOF Q-binning using grid method...\n')
        
        sw_status(0,1,fid);
        
        for ii = 1:prod(nQ)
            hklC = bsxfun(@plus,hkl,[dhkl{1}(ii); dhkl{1}(ii); dhkl{1}(ii)]);
            spec0 = obj.spinwave(hklC,'formfact',spec.formfact);
            spec0 = sw_egrid(spec0,'Evect',Evect,'component',spec.component);
            conv0 = conv0 + spec0.swConv;
            
            sw_status(ii/prod(nQ)*100,0,fid);
        end
        
    otherwise
        error('sw_tofres:WrongOption','Wrong method option!')
end

sw_status(100,2,fid);
obj.fileid(fid);
spec.swConv = conv0/prod(nQ);
fprintf0(fid,'Calculation finished.\n')

end