function d2ddat = corrphon(d2ddat, varargin)
% subtract incoherent phonons from magnetic data in d2d object
%
% NDBASE.CORRPHON(d2d,'option1',value1,...)
%
% Options:
%
% qphon         Q-range where the phonon background is to be fitted in
%               units of A^-1. Recommended range is above 4 A^-1, where the
%               magnetic signal is negligible due to the form factor. It is
%               also possible to given an upper value. For no upper limit,
%               use inf. Default value is [4 5].
%

inpForm.fname  = {'qphon'};
inpForm.defval = {[4 5]  };
inpForm.size   = {[1 2]  };

param = sw_readparam(inpForm, varargin{:});

if ~isa(d2ddat,'d2d')
	error('corrphon:WrongInput','The given data is not d2d type!')
end

% convert data, x-axis in A^-1 units
dat = ndbase.getdata(d2ddat,'binType','center','axUnit','A-1');

% data to fit
sigFit = dat.sig(dat.ax{1}>param.qphon(1) & dat.ax{1}<param.qphon(2),:);

end