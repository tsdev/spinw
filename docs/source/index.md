---
title: Documentation for SpinW
keywords: homepage
sidebar: sw_sidebar
permalink: index.html
summary:
mathjax: true
---

This is the **official documentation** of SpinW. The documentation contains the general description of the code (under the "Documentation" menu) and the function reference for all SpinW classes, packages and function.

SpinW is a Matlab library that can define spin Hamiltonians including general quadratic exchange interactions, single ion anistropy and external magnetic field, applying space group symmetries and solve the resulting equation of motion.

SpinW constist of a [spinw](spinw) class, that stores all parameters of the spin Hamiltonian as a set of properties and has several methods to prepare the Hamiltonian and calculate the spin-spin correlation function.

The [swfiles](swfiles) folder contains additional function that can be used to post-process the calculated spin-spin correlation function, e.g. calculating cross section or convoluting instrumental resolution.

The [swplot](swplot) package contains functions that can visualize the spin Hamiltonian as 3D plot including bonds, magnetic moments and exchange interactions.

The [swpref](swpref) package stores persistent preferences, [swsym](swsym) package manipulates space group operators and [swfunc](swfunc) package contains spectral peak shape functions.

**Feel free to ask questions & requests!**


{% include links.html %}