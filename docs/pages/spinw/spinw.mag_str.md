---
{title: spinw.mag_str property, link: spinw.mag_str, summary: stores the magnetic
    structure, keywords: sample, sidebar: sw_sidebar, permalink: spinw_mag_str.html,
  folder: spinw, mathjax: 'true'}

---
 
### Sub fields
 
`F`
: Complex magnetization (strictly speaking complex
  spin expectation value) for every spin in the magnetic
  cell, represented by a matrix with dimensions of $$[3\times
  n_{magext}\times n_k]$$,
  where `nMagExt = nMagAtom*prod(N_ext)` and $$n_k$$ is the number
  of the magnetic propagation vectors.
 
`k`
: Magnetic propagation vectors stored in a matrix with dimensions
  of $$[3\times n_k]$$.
 
`N_ext`
: Size of the magnetic supercell in lattice units, default value
  is `[1 1 1]` emaning that the magnetic cell is identical to the
  crystallographic cell. The $$[1\times 3]$$ vector extends the cell
  along the $$a$$, $$b$$ and $$c$$ axes.
 

