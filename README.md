# 2D Lid-Driven Cavity Flow Simulation

This project implements a 2D lid-driven cavity flow simulation using the finite difference method. The solution is based on the streamfunction-vorticity formulation of the incompressible Navier–Stokes equations.

## Tech Stack

- C Language  
  Used to implement the numerical algorithm and solve the flow field.

- Finite Difference Method

- Tecplot  
  Used to visualize streamlines, vorticity contours, velocity vectors and velocity profiles.

- File Formats  
  `.c` source code, `.dat` output files.

## Simulation Details

- Grid size: 128 × 128
- Reynolds numbers: 100 and 400
- Outputs:
  - Streamfunction
  - Vorticity
  - Velocity profiles (u and v)
  - Visualization-ready `.dat` files for Tecplot

## Validation

Results are validated by comparing velocity profiles with benchmark data from Ghia et al. (1982), showing good agreement.

## Output Files

- streamlines.dat  
- vorticity.dat  
- u_velocity.dat  
- v_velocity.dat  

These files can be opened in Tecplot for plotting.

## Reference

Ghia, U., Ghia, K. N., & Shin, C. T. (1982).  

## Author

Neeraj Kumar


