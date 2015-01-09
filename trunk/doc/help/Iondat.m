%% Ion.dat
% The ion.dat file contains information about different magnetic ions,
% where each row defines an ion. The meaning of the columns are
% the following:
%
% * *Name* of ion, string. These labels can be used in the
% sw.unit_cell.label field to assign different ions to crystallographic
% positions, predefined names follow the [M][element name][charge] notation
% similarly to FullProf.
% * *Magnetic form factor*, double. It is either 7 or 9 doubles the A, a,
% B, b, ... coefficients in the formula:
% LATEX\langle j_0(Q_s)\rangle = A\cdot exp(-a\cdot Q_s^2) + B\cdot exp(-b\cdot Q_s^2) + C* exp(-c\cdot Q_s^2)+D\cdot exp(-d\cdot Q_s^2) + EPATEX.
% * *Spin quantum number*, double. This number defines the spin value that
% is used as default in the sw.addatom mathod if 'S' is undefined.
%
% See also atom.dat, SW_MFF.