---
{title: sw_fstat, link: sw_fstat, summary: calculates termodynamical averages during
    an annealing simulation, keywords: sample, sidebar: sw_sidebar, permalink: sw_fstat.html,
  folder: swfiles, mathjax: 'true'}

---

### Syntax

`parout = sw_fstat(state, parin, t, e, m, next)`

### Description

It calculates statistical properties of different physical variables over
several sample state. Called by spinw.anneal.
 

### Input Arguments

`state`
: Defines the task of the function.
  1   Initialize the parOut structure.
  2   Store the parameters of the physical state.
  3   Calculate physical properties from the variable
      statistics.

`parIn`
: Same as parOut.

`T`
: Temperature of the system, vector [1,nT].

`E`
: Energy of the system, vector [1 nT].

`M`
: Magnetic moment of every atom, size: [spinDim,nMagExt*nT].

`nExt`
: Size of the magnetic cell, size: [3,1].

`kB`
: Boltmann constant, units of temperature.

### Output Arguments

parOut        Output parameter structure.
parOut.nStat  The number of evaluated states.
parOut.M      <M> summed over all magnetic moment, dimensions are
              [spinDim,nMagExt*nT].
parOut.M2     <M^2> summed over all magnetic moment, dimensions are
              [spinDim,nMagExt*nT].
parOut.E      <E> summed over all magnetic moment.
parOut.E2     <E^2> summed over all magnetic moment.
For the final execution, the following parameters are calculated:
parOut        Array of struct, size [1 nT].
parOut.avgM   Average components of the magnetisation over nStat runs,
              size: (3,nMagExt).
parOut.stdM   Standard deviation of the mgnetisation components over
              nStat runs, size: (3,nMagExt).
parOut.avgE   Average system energy per spin over nStat runs, scalar.
parOut.stdE   Standard deviation of the system energy per spin over
              nStat runs, scalar.
parOut.T      Final temperature of the sample.
parOut.Cp     Heat capacity of the sample: (<E^2>-<E>^2)/kB/T^2.
parOut.Chi    Magnetic susceptibility of the sample: (<M^2>-<M>^2)/kB/T.

### See Also

[spinw](spinw.html) \| [spinw.anneal](spinw_anneal.html)

