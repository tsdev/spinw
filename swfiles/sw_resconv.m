function M = sw_resconv(M,x,dx,func0)
% Convolute Gaussian with variable width along the first dimension of a matrix
%
% M = SW_RESCONV(M,x,dx,func)
%
% This convolution keeps the integrated intensity (sum(I*dx)) constant. It
% assumes the x vector contains the center points of the bins and the
% distances between the generated bin edges is calculated by interpolating
% from the distances between the given x bin center positions.
%
% Input:
%
% M     Arbitrary matrix with dimensions of (m1,m2).
% x     Column vector of coordinates along the first dimension of the
%       matrix.
% dx    Full width at half maximum (FWHM) value of the Gaussian as a
%       function of dx. Either a function handle with a header:
%           fwhm = dx(xVal)
%       or a vector of polynomial coefficients that produces the right
%       standard deviation. In this case in the function the following line
%       will be executed:
%           fwhm = polyval(dx,xVal)
%       The standard deviation of the Gaussian is calculated from the given
%       dx value using the following formula:
%           stdG = fwhmG/2/sqrt(2*log(2)) ~ fwhmG/2.35482
%       If a general resolution function is provided in the func argument,
%       it will be called as:
%           y = func(x,[1 x0 fwhm]);
%       In this case the fwhm can be a row vector and the meaning of the
%       different parameters will depend on func.
% func  Resolution function shape with header:
%           y = func(x,p)
%       Where x is a column vector, p is a row vector of parameters. The
%       meaning of the first 2 elements of the parameter vector are fixed.
%           p(1)    integral of the function.
%           p(2)    center of mass position of the function.
%       Optional, default value is @swfunc.gaussfwhm.
%
% Output:
%
% M     Matrix with same dimensions as the input storing the convoluted
%       data.
%
% See also SW_RES.
%

if nargin == 0
    help sw_resconv
    return
end

if nargin < 4
    func0 = @swfunc.gaussfwhm;
end

Mtemp = M * 0;

% calculate bin size from center bins (x)
bin = diff(x(:))';
bin = [bin(1) (bin(1:(end-1))+bin(2:end))/2 bin(end)];

for ii = 1:numel(x)
    % standard deviation of the energy resolution Gaussian
    if isa(dx,'function_handle')
        fwhm = dx(x(ii));
    else
        fwhm = polyval(dx,x(ii));
    end
    
    % Gaussian with intensity normalised to 1, centered on E(ii)
    %fG = 1/sqrt(2*pi)/stdG*exp(-((x-x(ii))/stdG).^2/2);
    % proper normalization should work also for unequal bins
    fG = func0(x,[1 x(ii) fwhm])*bin(ii);
    
    Mtemp = Mtemp + fG * M(ii,:);
    
end

M = Mtemp;

end