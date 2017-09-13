---
{title: sw_nvect( ), summary: determines the best normal vector for the set of vectors,
  keywords: sample, sidebar: sw_sidebar, permalink: sw_nvect.html, folder: swfiles,
  mathjax: 'true'}

---
determines the best normal vector for the set of vectors
 
[n collinear] = SW_NVECT(S, {epsilon})
 
The function can deal with complex vectors, separating the real and
complex parts as separate vectors.
 
S           Array of column vectors, dimensions are [3 N].
epsilon     Upper limit of the collinearity and the lower limit of
            coplanarity. Default value is 0.1. If epsilon = 1 the
            function returns n vector closest to the collinear direction,
            if epsilon = 2, n will be closest to the normal of the spin
            plane.
 
n           Vector parallel to the collinear spin direction or
            perpendicular to the plane of the spins in planar structures,
            dimensions are [1 3].
 
collinear   If true, the set of vectors are collinear, then n is parallel to the
            input vectors.
 
