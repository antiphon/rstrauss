# rstrauss: R-package for simulating and fitting Strauss point pattern model in 2- and 3-dimensions (2D, 3D)

Simulate and fits the stationary and isotropic Strauss point process in 2D or 3D cube (rectangular cuboid). 

Implements standard simulation algorithms: Metropolis-Hastings, birth-and-death and dominated CFTP.

Implements the approximate maximum logistic likelihood estimation method and direct numerical maximum likelihood estimation using approximations of the normalizing constant.

Partial support also for higher than 3D.

Changes:

- v1.3: Remove depenedency on vblogistic. Minor
  - 1: fstrauss log-likelihood correction was missing in sampling based fstrauss
  
  
  
todo:
  - add support for non-rectangular windows