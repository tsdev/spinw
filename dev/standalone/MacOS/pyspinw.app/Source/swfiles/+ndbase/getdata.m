function dat = getdata(obj, varargin)
% extract data from different type of objects
%
% dat = NDBASE.GETDATA(dndobj, 'option1',value1,...)
%
% Currently works with dnd objects.
%
%
% Input:
%
% dndobj        Object of dnd type.
%
% Options:
%
% binType       Type of the output bin:
%                   'center'    Center bin, default.
%                   'edge'      Edge bin (storing edge value of each bin).
% axUnit        Unit of the axis:
%                   'default'   Units as it is given in the dnd object,
%                               default.
%                   'A-1'       Units of A^-1.
% emptyval      Value for empty pixels. Default is nan.
%
% Output:
%
% dat           Structure with the following fields:
%                   sig     Storing the signal in a matrix with dimensions
%                           [nAx1 nAx2 ...].
%                   err     Standard deviation of the signal stored in a
%                           matrix with same dimensions.
%                   ax      Axes stored in a cell: {ax1 ax2 ...}. Each axn
%                           is a column vector with nAxn or nAxn+1 number
%                           of elements for center and edge bins
%                           respectively.
%

inpForm.fname  = {'bintype' 'axunit'  'emptyval'};
inpForm.defval = {'center'  'default' nan       };
inpForm.size   = {[1 -1]    [1 -2]    [1 1]     };

param = sw_readparam(inpForm, varargin{:});

switch param.bintype
    case {'center' 'centre'}
        cbin = true;
    case 'edge'
        cbin = false;
    otherwise
        error('getdata:WrongOption','The given binType option is invalid!')
end

switch param.axunit
    case 'default'
        convax = false;
    case {'A-1' 'a-1'}
        convax = true;
    otherwise
        error('getdata:WrongOption','The given axUnit option is invalid!')
end

nDat = numel(obj);

dat = struct('ax',cell(1,nDat),'sig',cell(1,nDat),'err',cell(1,nDat));

for ii = 1:numel(obj)
    objS = struct(obj(ii));
    
    % get signal and error (use
    dat(ii).sig = objS.s;
    dat(ii).sig(objS.npix==0) = param.emptyval;
    dat(ii).err = sqrt(objS.e);
    dat(ii).err(objS.npix==0) = param.emptyval;
    
    % get axes
    dat(ii).ax = objS.p;
    
    if cbin
        % convert to center bin
        for jj = 1:numel(dat(ii).ax)
            dat(ii).ax{jj} = (dat(ii).ax{jj}(1:(end-1))+dat(ii).ax{jj}(2:end))/2;
        end
    end
    if convax
        % convert to A-1 units
        for jj = 1:numel(dat(ii).ax)
            dat(ii).ax{jj} = dat(ii).ax{jj}*objS.ulen(objS.pax(jj));
        end
        
    end
end

end