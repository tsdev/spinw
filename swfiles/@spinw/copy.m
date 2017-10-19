function objC = copy(obj)
% clones spinw object
% 
% ### Syntax
% 
% `newObj = copy(obj)`
% 
% ### Description
% 
% `newObj = copy(obj)` clones a SpinW object with all of its internal data.
% The `newObj` will be independent of the original `obj`. Since the [spinw]
% is a handle class, this command should be used to duplicate an object
% instead of the `=` operator. Using the `=` operator does not create a new
% object, but only a pointer that points to the original object:
% ```
% obj2 = obj
% ```
% Changing `obj` after the above command will also change `obj2`.
%
% ### Examples
% 
% In this example $J_{1a}$ is a matrix with 1 in the diagonal, while
% $J_{1b}$ has 3.1415 in the diagonal. If `cryst` is changed, `cryst1` will
% be changed as well and viece versa, since they point to the
% same object. However `cryst2` is independent of `cryst`:
%
% ```
% >>cryst = spinw
% >>cryst.addmatrix('label','J1','value',3.1415)
% >>cryst1 = cryst
% >>cryst2 = cryst.copy
% >>cryst.addmatrix('label','J1','value',1)
% >>J1a = cryst1.matrix.mat>>
% >>J1b = cryst2.matrix.mat>>
% ```
%
% ### Input Arguments
% 
% `obj`
% : [spinw] object.
% 
% ### Output Arguments
% 
% `newObj`
% : New [spinw] object that contains all the data of `obj`.
% 
% ### See Also
% 
% [spinw] \| [spinw.struct]
%

objS = struct(obj);
objC = spinw(objS);

% copy the private properties
objC.propl  = event.proplistener.empty;
objC.sym    = obj.sym;
objC.symb   = obj.symb;
objC.ver    = obj.ver;

% add new listeners to the new object
if ~isempty(obj.cache.matom)
    % add listener to lattice and unit_cell fields
    objC.addlistenermulti(1);
end

if ~isempty(obj.cache.symop)
    % add listener to lattice, unit_cell and coupling fields
    objC.addlistenermulti(2);
end

end