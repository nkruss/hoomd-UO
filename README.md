Author: Noah Kruss

Group: Paulose Group (University of Oregon - Physics)

# Molecular Dynamics Simulation and Analysis Project for Spring-Mass Systems
# Possessing time-dynamical springs

## System Description

This system is a code repository for running and analyzing molecular dynamic
simulations, preformed using the Hoomd-Blue software package, of spring-mass
systems of various lattice configurations possessing springs whose stiffnesses
vary sinusoidally in time.  

## Software Dependencies

Hoomd-Blue (https://hoomd-blue.readthedocs.io/en/latest/tutorial/00-Introducing-HOOMD-blue/00-index.html)

## Repo Organization

* analysis_files - A directory containing files that contain analysis functions
  for completed simulations

* aural_metamaterial - A directory to store the initial conditions of particle
  positions and velocities for the aural meta-material project along with simulation
  output folders

* sandbox - A directory containing files that contain various analysis and
  simulation functions that were used for various in-progress testing

* simulation_files - A directory containing files for initializing various
  lattice configurations within the Hoomd-Blue framework, creating dynamic
  springs, and running the molecular dynamic simulations

* 1D_line_test.py - File that contains functions for running simulations and
preforming analysis on a 1D spring mass chain with one of the following conditions
    - Being forced from one end at a set frequency
    - Having a 1D standing wave initialized within the chain

* aural_metamaterial.py - Functions for running the aural meta-material
  simulations

* aural_paper_figures.py - File that contains functions for creating the final
  figures for the aural meta-material paper
