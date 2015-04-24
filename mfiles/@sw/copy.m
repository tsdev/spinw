function objC = copy(obj)
% clones sw object with all data
%
% newObj = COPY(obj)
%
% Use this command instead of the '=' sign if you want to
% create an independent duplicate of an sw class object.
%
% Input:
%
% obj       sw class object.
%
% Output:
%
% newObj    New sw class object that contains all the data of
%           obj.
%
% Example:
%
% cryst = sw;
% cryst.addmatrix('label','J1','value',3.1415)
%
% cryst1 = cryst;
% cryst2 = cryst.copy;
%
% cryst.addmatrix('label','J1','value',1)
% J1a = cryst1.matrix.mat;
% J1b = cryst2.matrix.mat;
%
% Where J1a will be a matrix with 1 in the diagonal, while J1b
% has 3.1415 in the diagonal. If cryst is changed, cryst1 will
% be changed as well and viece versa, since they point to the
% same object. However cryst2 is independent of cryst.
%
% See also SW, SW.STRUCT.
%

objS = struct(obj);
objC = sw(objS);

% copy the private properties
objC.sym    = obj.sym;
objC.symb   = obj.symb;
objC.fid    = obj.fid;
objC.Elabel = obj.Elabel;
objC.Qlabel = obj.Qlabel;
objC.Rlabel = obj.Rlabel;
objC.Blabel = obj.Blabel;
objC.Tlabel = obj.Tlabel;

end