# N-Body Chaos Simulation

A computational physics project developed to investigate chaotic behavior in self-gravitating N-body systems.

## Overview

This project simulates the evolution of a gravitational system composed of multiple interacting particles using a direct N-body approach and Verlet integration.

The main objective is to study how extremely small differences in the initial conditions can lead to significantly different trajectories over time, illustrating the phenomenon of deterministic chaos.

The simulation also tracks several physical quantities throughout the evolution of the system, allowing qualitative analysis of its dynamical and thermodynamical behavior.

## Physical Concepts

* Newtonian Gravitation
* N-Body Dynamics
* Deterministic Chaos
* Virial Equilibrium
* Energy Conservation
* Statistical Temperature
* Half-Mass Radius

## Numerical Methods

The equations of motion are integrated using the Verlet algorithm, a widely used method in computational physics due to its numerical stability and good conservation properties.

The gravitational force is softened through a softening parameter in order to avoid numerical divergences during close encounters.

## Diagnostics

The code computes:

* Total Energy
* Kinetic Energy
* Potential Energy
* Virial Parameter
* Statistical Temperature
* Half-Mass Radius
* Trajectory Divergence Between Two Nearly Identical Systems

The divergence between two simulations differing only by a tiny perturbation in the initial conditions is used as a simple indicator of chaotic behavior.

## Project Structure

src/

* N_corpos.cpp

data/

* caos.csv

## Author

Eduardo Alebrandt Pires

Physics Engineering Student — UFRGS