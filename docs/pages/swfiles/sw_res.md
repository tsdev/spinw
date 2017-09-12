---
{title: sw_res( ), summary: reads a tabulated energy resolution from a file and fits
    with polynomial, keywords: sample, sidebar: sw_sidebar, permalink: swfiles_sw_res.html,
  folder: swfiles, mathjax: 'true'}

---
reads a tabulated energy resolution from a file and fits with polynomial
 
p = SW_RES(fid,polDeg,{plot})
 
The file contains the FWHM energy resolution values as a function of
energy transfer in two columns, first the energy transfer values
(positive is energy loss), the second is the FWHM of the Gaussian
resolution at the given energy transfer value.
 
Input:
 
fid           String, path to the resolution file or a matrix with the
              same format as the data file.
polDeg        Degree of the fitted polynomial to the instrumental
              resolution data. Default is 5.
plot          If true the resolution will be plotted, optional, default
              is true.
 
Output:
 
p             Returns the coefficients for a polynomial p(x) of degree n
              that is a best fit (in a least-squares sense) for the resolution data
              in y. The coefficients in p are in descending powers, and
              the length of p is n+1.
              p(x)=p_1*x^n+p_2*x^(n-1)+...+p_n*x+p_(n+1).
 
Example:
 
To calculate the resolution at an arbitrary energy use:
 
pRes = sw_res(file,5,false);
Evec = linspace(0,100,501);
Eres = polyval(polyRes,Evec);
plot(Evec,Eres);
 
 
See also POLYFIT, SW_INSTRUMENT.
 
