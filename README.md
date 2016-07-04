[![Build Status](http://ci.sagrid.ac.za/buildStatus/icon?job=Fermionic-Molecular-Dynamics-deploy)](http://ci.sagrid.ac.za/job/Fermionic-Molecular-Dynamics-deploy)

# FMD codes (from  T. Neff)

This is a repository of codes NOT CREATED BY ME, but by my advisor, Thomas Neff. They are for FMD calculations.

# Building

## Quickstart

For the impatient:

  1. Edit `Makefile.inc` based on one of the tempates provided, to point to your local compiler and dependencies
  2. Run `make`
  3. SCIENCE !


## Dependencies

The application requires :

  * GCC, gfortran
  * ncurses
  * LAPACK (libblas, liblapack)
  * OpenMPI

# Testing and Delivery

This repository is automatically integrated with [Project CODE-RADE](http://www.africa-grid.org/CODE-RADE). Build, test and deploy scripts are in the [CODE-RADE](CODE-RADE) folder.

  * `build.sh` - this script runs make and builds the static libraries as well as the executable files.
  * `check-build.sh` - this script checks whether the build is successful by running tests
  * `deploy.sh` - this script cleans the build, reconfigures it for deployment (to `/cvmfs`) and then rebuilds with the new target.
