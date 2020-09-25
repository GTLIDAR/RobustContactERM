# RobustContactERM
Repository for Robust Trajectory Optimization with Stochastic Complementarity

**[Release Note Setpember 2020]** 
This work has been submitted to RA-L and to arXiv *link pending*

## Overview
This repository contains code for robust contact-implicit trajectory optimization in the presence of uncertainty in the terrain model, as described in the paper *Robust Trajectory Optimization over Uncertain Terrain with Stochastic Complementarity*. Our implementation uses Expected Residual Minimization, a smoothing method for solving complementarity problems with uncertain data, and is tested in three examples including a footed hopping robot and a benchmarking experiment on sliding a block over terrain with uncertain friction characteristics.

* Authors: Luke Drnach and Ye Zhao
* Affiliation: [The LIDAR Lab](http://lab-idar.gatech.edu/), Georgia Institute of Technology

This code was tested in *MATLAB 2017a* with *Ubuntu 16.04* and a MATlAB version of [*Drake*](https://drake.mit.edu)

