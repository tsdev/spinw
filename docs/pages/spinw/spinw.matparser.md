---
{title: spinw.matparser method, link: spinw.matparser, summary: assigns new values
    to existing matrices, keywords: sample, sidebar: sw_sidebar, permalink: spinw_matparser.html,
  folder: spinw, mathjax: 'true'}

---
 
MATPARSER(obj, 'option1', value1 ...)
 
The function modifies the spinw.matrix.mat matrix, it assigns new values
from a given parmeter vector.
 
Input:
 
obj           Input structure, spinw class object.
 
Options:
 
param         Input vector P with nPar elements that contains the
              new values to be assignd to elements of spinw.matrix.mat
              matrix.
mat           Identifies which matrices to be changed according to their
              label or index. To select matrices with given labels use a
              cell of strings with nPar elements, for example
              M = {'J1','J2'}. This will change the diagonal elements of
              matrices J1 and J2 to a given value that is provided in the
              'param' option. Alternatively the index of the matrices can
              be given in a vector, such as [1 2] (index runs according
              to the order of the previous creation of the matrices using
              spinw.addmatrix).
              To assign parameter value only to a selected element of a
              3x3 matrix, a bracket notation can be used in any string,
              such as 'D(3,3)', in this case only the (3,3) element of
              the 3x3 matrix of 'D' will be modified, the other elements
              will be unchanged. To modify multiple elements of a matrix
              at once, use the option 'selector'.
selector      Matrix with dimensions of [3 3 nPar]. Each S(:,:,ii)
              submatrix can contain +/-1 and 0. Where S(:,:,ii) contains
              ones, the corresponding matrix elements of
              spinw.matrix.mat(:,:,M(ii)) will be changed to the value
              P(ii)*S(:,:,ii) where P(ii) is the corresponding parameter
              value. For example do assign a Dzyaloshinskii-Moriya vector
              to the 'DM' matrix, the following input would be
              sufficient:
              P = [0.2 0.35 3.14]
              M = {'DM' 'DM' 'DM'}
              S = cat(3,[0 0 0;0 0 1;0 -1 0],[0 0 -1;0 0 0;1 0 0],[0 1 0;-1 0 0;0 0 0])
              spinw.matparser('param',P,'mat',M,'selector',S)
init          Initialize the matrices of spinw.matrix.mat with zeros for all
              selected labels before assigning parameter values. Default
              is false.
 
Output:
 
The spinw object will contain the modified matrix.mat field.
 
See also SPINW, SPINW.HORACE, SPINW.ADDMATRIX.
 

