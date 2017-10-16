---
title: Introduction
keywords: docs
sidebar: sw_sidebar
permalink: crystal.html
summary: Quick introduction to SpinW
folder: documentation
mathjax: 'true'
---


## Requirements and Installation

### What you need

SpinW can be used in two different ways:
* The best is to run it from a recent version of [Matlab](https://www.mathworks.com/products/matlab.html) (newer than R2014a), however this requires a valid license.
* The alternative way is to run pySpinW from [Python](https://www.python.org) which does not require Matlab license and it can also be used for distributed computing.

The symbolic calculation mode also needs the [Symbolic Math Toolbox](https://www.mathworks.com/products/symbolic.html) installed.

### Installation

The latest version of SpinW can be downloaded from [https://github.com/tsdev/spinw/releases](https://github.com/tsdev/spinw/releases. Steps to install:
* extract the `.zip` file into any local folder
* add the folder to the Matlab search path using `addpath(genpath(MYSPINWFOLDER))`
* recommended to copy the above command into the `startup.m` file

### Updating

Since SpinW is under active development, there are regular updates (new features, bug fixes). To update your local copy of SpinW do the following:
* execute the [sw_update](sw_update) function in the Matlab Command Window, this will check for the latest available version and downloads it into a new subfolder under /spinw/ and adds the new files to the Matlab search path
* to make the new path permanent, edit the `startup.m` file as written above
* execute the `clear classes` command in the Command Window to clear the Matlab cache of user defined classes


## Hello World!

It is possible to create, solve and plot a spin wave model in a single line of code using SpinW! For example the following line will colve the spin wave dispersion of the triangular lattice antiferromagnet along the $(h,h,0)$ direction in reciprocal space:

```
>>sw_plotspec(spinwave(sw_model('triAF',1),{[0 0 0] [1 1 0]}))
>>snapnow
```



## Get help

## Handle Class

## Meta Rules

## SpinW Cheatsheet

{% include links.html %}
