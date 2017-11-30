function [polyRes, yout] = sw_res(fid,polDeg,toplot,varargin)
% fits energy resolution with polynomial
% 
% ### Syntax
% 
% `p = sw_res(source,poldeg)`
% 
% `p = sw_res(source,poldeg,true)`
%
% ### Description
% 
% `p = sw_res(fid,poldeg)` reads tabulated resolution data from the
% `source` file which contains the FWHM energy resolution values as a
% function of energy transfer in two columns. First  column is the energy
% transfer values (positive is energy loss), while the second is the FWHM
% of the Gaussian resolution at the given energy transfer.
% 
% `p = sw_res(fid,poldeg,plot)` the polynomial fit will be shown in a
% figure if `plot` is true.
%
% ### Examples
% 
% This example shows how to fit a tabulated resolution data (MERLIN energy
% resolution for $E_i=50$ meV and 300 Hz chopper frequency). Using the
% fitted polynomial, the energy resolution can be calculated at an
% arbitrary energy transfer value.
%
% ```
% >>resDat = [0 2.31;10 1.80;20 1.37;30 1.02;40 0.78;49 0.67]>>
% >>polyRes = sw_res(resDat,3)
% >>snapnow
% >>EN = 13
% >>dE = polyval(polyRes,EN)>>
% ```
% 
% ### Input Arguments
% 
% `source`
% : String, path to the resolution file or a matrix with two columns, where
%   the first column gives the energy transfer value and second column
%   gives the resolution FWHM.
% 
% `polDeg`
% : Degree of the fitted polynomial, default value is 5.
% 
% `plot`
% : If `true` the resolution will be plotted, default value
%   is `true`.
% 
% ### Output Arguments
% 
% `p`
% : The coefficients for a polynomial $p(x)$ of degree $n$
%   that is a best fit (in a least-squares sense) for the resolution data.
%   The coefficients in $p$  are in descending powers, and
%   the length of $p$ is $n+1$:
%
%   $p(x)=p_1\cdot x^n+p_2\cdot x^{n-1}+...+p_n\cdot x+p_{n+1}$
% 
% ### See Also
% 
% [polyfit] \| [sw_instrument]
%
% *[FWHM]: Full Width at Half Maximum
%

if nargin == 0
    swhelp sw_res
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