# Robust_tensor_recovery_for_traffic_events
Yue Hu, Jun 2019

## Overview
This repository contains the source code developed for robust tensor recovery for traffic events detection. This results are reported in "Robust Tensor Recovery with Fiber Outliers for TrafficEvents" by Yue Hu and Daniel Work, which was submitted to ACM Transactions on Knowledge Discovery from Data (TKDD).

## Structures
- `/code/` The source code folder.
  - `inexact_alm_rpca21.m` and `inexact_alm_rpca3D.m` deals with full observation case, and sovle robust tensor PCA under fiber-sparse corruption and element-sparse corrution, respectively.
  - `inexact_alm_rmc21.m` and `inexact_alm_rmc3D.m` deals with partial observation case, and solve tensor robust complepetion under fiber-sparse corruption and element-sparse corrution, respectively.
  -  `inexact_alm_rmc2D_21.m` solves matrix robust completion under column-wise corruption.
  - `/test/` Test code folder, including simulation tests and traffic data test.
  - `/PROPACK/` Prerequisit packages, including code for efficient PCA and tensor manipulation package - `Tensor Toolbox for MATLAB`.
- `/SimulationData/` Contains results of simulation tests,
- `/Figure/` Contains figures for traffic tensor arrangement. Simulation result figures are also saved here. 
- `/TrafficData/` The traffic dataset folder, which contains traffic observation data arranged in tensor format.


## Contact
+ Author: Yue Hu, Institute for Software Integrated Systems, Vanderbilt University
+ Email: yue.hu(at)vanderbilt.edu
