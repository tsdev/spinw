---
{title: sw_parstr( ), link: sw_parstr, summary: parses input string, keywords: sample,
  sidebar: sw_sidebar, permalink: sw_parstr.html, folder: swfiles, mathjax: 'true'}

---

### Syntax

` `

### Description

sections in strIn.
 

### Input Arguments

% `'Sperp'`
:   convolutes the magnetic neutron scattering intensity

% `(<Sperp`
: * Sperp> expectation value). Default.

% `'Sab'`
:   convolutes the selected components of the spin-spin correlation

% `function.`
:n. Letter a and b can be 'x', 'y' or 'z'. For example:

% `'Sxx'`
:ill convolute the xx component of the correlation

% `function`
:n with the dispersion. xyz is the standard coordinate

% `system,`
: see online documentation of spinw.

% `'Mab'`
:   convolutes the selected components of the spin-spin

% `correlation`
:tion function. Letter a and b can be 'x', 'y' or 'z'.

% `For`
:mple: 'Sxx' will convolute the xx component of the

% `correlation`
:tion function with the dispersion. The xyz coordinates

% `are`
:the Blume-Maleev coordinate system, see below.

% `'Pab'`
:   convolutes the selected element of the polarisation

% `matrix.`
: Letter a and b can be 'x', 'y' or 'z'. For example:

% `'Pyy'`
:ill convolute the yy component of the polarisation

% `matrix`
:with the dispersion. The xyz coordinates are in the

% `Blume-Maleev`
:aleev coordinate system, see below.

% `'Pa'`
:   convolutes the intensity of the simulated polarised

% `neutron`
: scattering, with inciden polarisation of Pa. Letter a

% `can`
:'x', 'y' or 'z'. For example: 'Py' will convolute the

% `scattering`
:ing intensity simulated for incident polarisation Pi ||

% `y.`
:xyz coordinates are in the Blume-Maleev coordinate

% `system,`
: see below.

% `Any`
:ear combination of the above are allowed, for example: 'Sxx+2*Syy'

% `convolutes`
:tes the linear combination of the xx component of the spin-spin

% `correlation`
:tion function and the yy component.

% `The`
:me-Maleev coordinate system is a cartesian coordinate system

% `with`
:BM, yBM and zBM) basis vectors as follows:

% `xBM`
:parallel to the momentum transfer Q,

% `yBM`
:perpendicular to xBM in the scattering plane,

% `zBM`
:perpendicular to the scattering plane.

### Output Arguments

parsed is struct type, contains the following fields:
type      Cell contains as many elements as many in the sum. Each element
is a vector as follows:
type{idx}(1)    Index of type of cross section:
                1   Sperp,
                2   Sab,
                3   Mab,
                4   Pab,
                5   Pa.
type{idx}(2:3)  Index of the component:
                1   x,
                2   y,
                3   z.
preFact   Vector contains the values of the prefactors in the sum.
string    Original input string.
Test it with:
<a href="matlab:parsed = sw_parstr('Sxx + Syy')">parsed = sw_parstr('Sxx + Syy')</a>

### See Also

[sw_egrid](sw_egrid.html) and [spinw.fitspec](spinw_fitspec.html)

