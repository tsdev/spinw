function sFact = sw_intsf(sFact, varargin)
% integrates the structure factor along given Q directions
%
% sFact = SW_INTSF(sFact, 'Option1', Value1, ...) 
%
% Options:
%
% axis      The index of the axis, along which the data is summed within
%           the given range. {1, 2, 3} is for {a*, b*, c*} respectively.
%           Zero is to powder average the data. Default is 1. More than one
%           axis can be given.
% range     Data range in r.l.u. along selected dimension to integrate,
%           default is the full data range.
% type      Type of data to be used for the integration.
%               'sf'    integrate the calculated structure factor. (default)
%               'perp'  integrate the perpendicular component of the structure
%                       factor to Q.
%
% Output:
%           sFact contains the following new fields:
%           int     The integrated data in a matrix.
%           hklint  A cell containing the three axis vectors in r.l.u. and
%
% Example:
%   sFact = cryst.structfact;
%   sFact = sw_intsf(sFact,'axis',[1 2],'range',[0 1; 0 2]);
%
%   The above code will integrate the structure factor along a* between 0
%   and 1, along b* between 0 and 2.
%
% See also SPINW, SPINW.STRUCTFACT, SW_PLOTSF.
%

if nargin == 0
    help sw_intsf
    return
end

inpForm.fname  = {'axis' 'range' 'type' 'dQ'  };
inpForm.defval = {1      []      'sf'    []   };
inpForm.size   = {[1 -1] [-1 2]  [1 -2]  [1 1]};
inpForm.soft   = {0      1       0       1    };

param = sw_readparam(inpForm, varargin{:});

switch param.type
    case 'sf'
        F2 = sFact.F2;
    case 'perp'
        F2 = sFact.perp;
    otherwise
        warning('sw_plotsf:WrongParam','param.type is wrong, using ''sf''!');
        F2 = sFact.F2;
end


% sum data
if param.axis == 0
    % Powder average the data
    % create hkl grid
    [hh, kk, ll] = ndgrid(sFact.hkl{:});
    % calculate Q abs in inverse Angstrom
    %[h k l] * 2*pi*inv(basisvector)
    hkl  = permute(cat(4,hh,kk,ll),[5 4 1 2 3]);
    hklA = 2*pi*mmat(hkl,inv(sFact.obj.basisvector));
    QA   = permute(sqrt(sum(hklA.^2,2)),[3 4 5 1 2]);
    % determine the maximum Q, where the sphere is properly covered
    Qmax = min([QA(1,1,end),QA(1,end,1),QA(end,1,1)]);
    % create Q vector with the given steps
    if isempty(param.dQ)
        param.dQ = Qmax/100;
    end
    % Q bin
    epsilon = 1e-5;
    Qbin = [0:param.dQ:Qmax Qmax+epsilon];
    nQ   = numel(Qbin);
    [~, idxQ] = min(abs(repmat(QA(:),[1 nQ])-repmat(Qbin,[numel(QA) 1]))'); %#ok<UDIM>
    sFact.int = accumarray(idxQ',F2(:),[nQ 1]);
    N   = accumarray(idxQ',F2(:)*0+1,[nQ 1]);
    N(N==0) = 1;
    sFact.int = sFact.int./N;
    sFact.int = sFact.int(1:end-1);
    sFact.hklAint = Qbin(1:end-1);
    
    
else
    sIdx = cell(1,3);
    idx = 1;
    for ii = 1:3
        if (any(param.axis==ii)) && (idx<=size(param.range,1))
            % select axis vector
            sAxis = sFact.hkl{param.axis};
            sIdx{ii} = (sAxis>=param.range(idx,1)) & (sAxis<=param.range(idx,2));
            idx  = idx + 1;
        else
            sIdx{ii} = true(1,size(F2,ii));
        end
        
        sFact.hklint{ii} = sFact.hkl{ii}(sIdx{ii});
    end
    sFact.int = F2(sIdx{:});
    for ii = 1:numel(param.axis)
        sFact.int = sum(sFact.int,param.axis(ii));
    end
end

sFact.axis  = param.axis;
sFact.range = param.range;

end