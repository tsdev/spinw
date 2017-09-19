---
{title: sw_mattype( ), link: sw_mattype, summary: determines the type of square input
    matrix, keywords: sample, sidebar: sw_sidebar, permalink: sw_mattype.html, folder: swfiles,
  mathjax: 'true'}

---

### Syntax

` `

### Description

vector with dimensions of [1 N].
 

### Input Arguments

% `mat`
: Matrix with dimensions of [3 3 N].

% `epsilon`
: Error bar on small matrix element. Defult is 1e-5.

% ``
:

### Output Arguments

type      1   Heisenberg exchange
2   Anisotropic exchange
3   DM interaction
4   General matrix
Also works on symbolic matrices, but keep all symbols real for consistent
result!

