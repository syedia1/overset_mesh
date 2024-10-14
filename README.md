Implementing overset methodology framework for numerically soving equations.

The method has been demonstrated by solving 2D heat conduction and 2D Navier Stokes examples on structured grids.

THe modules implemented are:

1. Overlapping grids hole cutting - using geometric intersection algorithms (implemented using Alternating Digital Trees)
2. Interpolation algorithm

The repository work is sequentially ordered as:
1. Strucutred grid 2D heat conduction solver
2. One-way coupling between overset grids (2D heat conduction, one way information transfer from donor mesh to reciever(finer) mesh)
3. Two-way coupling between overset grids (2D heat conduction)
4. Two-way coupling between overset grids (Navier Stokes)
