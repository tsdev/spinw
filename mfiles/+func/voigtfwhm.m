function y = voigtfwhm(x,p)
% normalized function that calculates the voigt profile using FWMH values
%
% y = FITFUN.VOIGTFWHM(x,p)
%
% The integral of the function is normalized assumind dx = 1.
%
% Input:
%
% x     Input coordinates where the function will be calculated.
% p     Parameters:
%           A = p(1) integral of the signal assumin dx = 1 (for different
%           dx multiply the amplitude with dx to keep the integral
%           constant).
%           x0 = p(2) peak center positions.
%           wG = p(3) FWHM of the Gaussian component.
%           wL = p(4) FWHM of the Lorentzian component.
%
% Conversion between different width:
% gamma parameter of the Lorentzian
% gamma = wL/2
% Standard deviation of the Gaussian:
% sigma = wG/sqrt(8*ln(2))
%
% See also FUNC.GAUSS, FUNC.GAUSSFWHM.
%

% Origian code of:
% 27-December-2013 N. Cherkasov
% Comments and questions to: n.b.cherkasov@gmail.com

A  = p(1);
x0 = p(2);
wG = p(3)/2;
wL = p(4)/2;

% converting to dimensionless coordinates
x = sqrt(log(2)).*(x-x0)./(wG);
y = sqrt(log(2)).*(wL/wG);

w = complexErrorFunction(x,y);
y = A*sqrt(log(2)/pi)/wG.*real(w);

end

function w = complexErrorFunction(x,y)
% complexErrorFunction  Calculation of complex error function using dimentionless coordinates
%
% [w] = complexErrorFunction(x,y)   Computes the complex error function
%   using the algorithm developed by Dr. F. Schreier and kindly presented
%   in Fortran. The function was rewriten to MATLAB by Dr. N. Cherkasov
%   For more details on algorithm see the publication:
%   F. Schreier: Optimized Implementations of Rational Approximations for the Voigt ane Complex Error Function.
%   J. Quant. Spectrosc. & Radiat. Transfer, 112(6), 1010?1025, doi 10.1016/j.jqsrt.2010.12.010, 2011.
%
%   Briefly, the algorithm is compiled from two:
%       for    large x+y     J  Humlicek, JQSRT 27, 437, 1982
%       for small x+y:    J.A.C. Weideman,  SIAM J. Numer. Anal. 31 (1994) pp. 1497-1518,  equation (38.I) and table I
%
% INPUT ARGUMENTS are dimentioneless coordinates x and y
%   x - array 1*N, and y - single variable
%
% OUTPUT
%   w - complex array 1*N
%
% The function was used for the deconvolution of IR spectra
% see the publication
%
% 27-December-2013 N. Cherkasov
% Comments and questions to: n.b.cherkasov@gmail.com

%   "Weideman" constants
% n=24;
l = 4.1195342878142354; % l=sqrt(n/sqrt(2.))  ! L = 2**(-1/4) * N**(1/2)
a = [-1.5137461654527820e-10,  4.9048215867870488e-09,  1.3310461806370372e-09, -3.0082822811202271e-08, ...
    -1.9122258522976932e-08,  1.8738343486619108e-07,  2.5682641346701115e-07, -1.0856475790698251e-06, ...
    -3.0388931839840047e-06,  4.1394617248575527e-06,  3.0471066083243790e-05,  2.4331415462641969e-05, ...
    -2.0748431511424456e-04, -7.8166429956142650e-04, -4.9364269012806686e-04,  6.2150063629501763e-03, ...
    3.3723366855316413e-02,  1.0838723484566792e-01,  2.6549639598807689e-01,  5.3611395357291292e-01, ...
    9.2570871385886788e-01,  1.3948196733791203e+00,  1.8562864992055408e+00,  2.1978589365315417e+00];
%   humlicek prbFct region I bounds
s15 = 15;

% left wing -- center
x12 = y - s15;
% 15-y   center -- right wing
x21 = -x12;

if y>s15 || x(1)>x21 || x(end)<x12
    % all points are in Humlicek's region I
    t = y - x*1i;
    w = (t/sqrt(pi))./(1/2 + t.*t);
else
    s  = abs(x) + y;
    ds = s>s15;
    w  = zeros(size(x));
    
    % s>s15
    t = y - x(ds)*1i;
    w(ds) = (t/sqrt(pi))./(1/2 + t.*t);
    % s<s15
    recLmZ  = 1./(l+y-x(~ds)*1i);
    t       = (l-y+x(~ds)*1i).*recLmZ;
    w(~ds)  = recLmZ.*(1/sqrt(pi) + 2*recLmZ.*polyval(a,t));
    
end

end


