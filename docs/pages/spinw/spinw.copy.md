---
{title: spinw.copy method, link: spinw.copy, summary: clones spinw object, keywords: sample,
  sidebar: sw_sidebar, permalink: spinw_copy, folder: spinw, mathjax: 'true'}

---
  
### Syntax
  
`newObj = copy(obj)`
  
### Description
  
`newObj = copy(obj)` clones a SpinW object with all of its internal data.
The `newObj` will be independent of the original `obj`. Since the [spinw](spinw)
is a handle class, this command should be used to duplicate an object
instead of the `=` operator. Using the `=` operator does not create a new
object, but only a pointer that points to the original object:
```matlab
obj2 = obj
```
Changing `obj` after the above command will also change `obj2`.
 
### Examples
  
In this example $$J_{1a}$$ is a matrix with 1 in the diagonal, while
$$J_{1b}$$ has 3.1415 in the diagonal. If `cryst` is changed, `cryst1` will
be changed as well and viece versa, since they point to the
same object. However `cryst2` is independent of `cryst`:
 
```matlab
cryst = spinw
cryst.addmatrix('label','J1','value',3.1415)
cryst1 = cryst
cryst2 = cryst.copy
cryst.addmatrix('label','J1','value',1)
J1a = cryst1.matrix.mat
```
*Output*
```
J1a =
     1     0     0
     0     1     0
     0     0     1
```
 
```matlab
J1b = cryst2.matrix.mat
```
*Output*
```
J1b =
    3.1415         0         0
         0    3.1415         0
         0         0    3.1415
```
 
 
### Input Arguments
  
`obj`
: [spinw](spinw) object.
  
### Output Arguments
  
`newObj`
: New [spinw](spinw) object that contains all the data of `obj`.
  
### See Also
  
[spinw](spinw) \| [spinw.struct](spinw_struct)
 

{% include links.html %}
