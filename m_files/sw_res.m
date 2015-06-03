function [polyRes, yout] = sw_res(fid,polDeg,toplot,varargin)
% reads a tabulated energy resolution from a file and fits with polynomial
%
% p = SW_RES(fid,polDeg,{plot})
%
% The file contains the FWHM energy resolution values as a function of
% energy transfer in two columns, first the energy transfer values
% (positive is energy loss), the second is the FWHM of the Gaussian
% resolution at the given energy transfer value.
%
% Input:
%
% fid           String, path to the resolution file or a matrix with the
%               same format as the data file.
% polDeg        Degree of the fitted polynomial to the instrumental
%               resolution data. Default is 5.
% plot          If true the resolution will be plotted, optional, default
%               is true.
%
% Output:
%
% p             Returns the coefficients for a polynomial p(x) of degree n
%               that is a best fit (in a least-squares sense) for the resolution data
%               in y. The coefficients in p are in descending powers, and
%               the length of p is n+1.
%               p(x)=p_1*x^n+p_2*x^(n-1)+...+p_n*x+p_(n+1).
%
% See also POLYFIT, SW_INSTRUMENT.
%

if nargin == 0
    help sw_res;
    return
end

% default polynom degree
if nargin == 1
    polDeg = 5;
end

% plot the fit by default
if nargin < 3
    toplot = true;
end

% load file and read values
if ischar(fid)
    res = importdata(fid,' ');
    res = res(:,[1 2])';
else
    res = fid';
end

xres = res(1,:);
yres = res(2,:);

% fit polynom to instrument energy resolution
polyRes = polyfit(xres,yres,polDeg);
xnew = linspace(min(xres),max(xres),500);
ynew = polyval(polyRes,xnew);

if nargin > 3
    yout = polyval(polyRes,varargin{1});
end

if toplot
    fig0 = gcf;
    figure(fig0);
    plot(xres,yres,'o-');
    hold all
    plot(xnew,ynew,'r-');
    xlabel('Energy Transfer (meV)');
    ylabel('FWHM energy resolution (meV)');
    title('Polynomial fit the instrumental energy resolution');
    legend('Tabulated resolution',sprintf('Fitted polynomial (degree %d)',polDeg));
    figure(fig0);
end

end