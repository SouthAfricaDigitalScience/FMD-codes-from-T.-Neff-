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
