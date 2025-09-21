# 2D Lid-Driven Cavity Flow Simulation

This project is a C program that simulates 2D lid-driven cavity flow using the finite difference method. It solves the incompressible Navier–Stokes equations using the streamfunction-vorticity formulation.

## Features

- Simulates 2D lid-driven cavity flow
- Finite difference method with central differencing
- Grid size: 128 × 128
- Reynolds numbers: 100 and 400
- Generates:
  - Streamlines
  - Velocity vectors
  - Vorticity contours
  - Velocity profiles

## Validation

Velocity profiles from the simulation are compared with benchmark results from Ghia et al. (1982). The results show good agreement, confirming the accuracy of the method.

